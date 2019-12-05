{-# LANGUAGE TypeApplications #-}
{-# OPTIONS_GHC -fno-warn-type-defaults #-}

module ScalarWave where

import Control.Exception (assert)
import Control.Monad.Loops
import qualified Data.Vector.Unboxed as U
import Data.VectorSpace
import Dual
import Equations
import Foreign
import qualified Numeric.LinearAlgebra as L
import SSSFN
import System.CPUTime

default (Int)

--------------------------------------------------------------------------------

bndw :: Floating a => Coord a -> a
bndw (Coord u v) = swamp * cos (pi / 4 * swk * t) * cos (pi / 4 * swk * r)
  where
    t = v + u
    r = v - u
    swk = 1
    swamp = 0.1

guess1 :: Floating a => Coord a -> State a
-- guess1 c = State 1 0 (bndw c)
guess1 c = State 1 0 0

guess :: (L.Container L.Vector a, Floating a) => Grid (State a)
guess = gsample (\(u, v) -> guess1 (Coord u v))

--------------------------------------------------------------------------------

eqnInt1 :: Fractional a => Eqn -> Coord a -> State a -> State a -> State a -> State a -> State a -> State a -> a
eqnInt1 eqn c s su sv suv suu svv =
  case eqn of
    EqnG01 -> eqnG01 c s su sv suv
    EqnG22 -> eqnG22 c s su sv suv
    EqnSW2 -> eqnSW2 c s su sv suv
    EqnG00 -> eqnG00 c s su suu
    EqnG11 -> eqnG11 c s sv svv

eqnAxis1 :: Fractional a => Eqn -> Coord a -> State a -> State a -> State a -> State a -> State a -> State a -> a
eqnAxis1 eqn c@(Coord u v) s@(State g q w) su@(State gu qu wu) sv@(State gv qv wv) suv@(State guv quv wuv) suu@(State guu quu wuu) svv@(State gvv qvv wvv) =
  case eqn of
    EqnG01 -> eqnAxisG01 c s su sv suv
    EqnG22 -> eqnAxisG22 c s su sv suv
    EqnSW2 -> eqnAxisSW2 c s su sv suv
    EqnG00 -> eqnAxisG00 c s su sv suu
    EqnG11 -> eqnAxisG11 c s su sv svv

eqnOrig1 :: Floating a => Eqn -> Coord a -> State a -> State a -> State a -> State a -> State a -> State a -> a
eqnOrig1 eqn c s@(State g q w) su sv suv suu svv =
  case eqn of
    EqnG01 -> g - 1
    EqnG22 -> q - 0
    EqnSW2 -> w - bndw c
    EqnG00 -> 0
    EqnG11 -> 0

-- i == 0, u == -1
eqnBndU1 :: (Eq a, Floating a) => Eqn -> Coord a -> State a -> State a -> State a -> State a -> State a -> State a -> a
eqnBndU1 eqn c@(Coord u v) s@(State g q w) _ sv _ _ svv =
  assert (u == -1) $
    case eqn of
      EqnG01 -> eqnG11 c s sv svv
      EqnG22 -> q - 0
      EqnSW2 -> w - bndw c
      EqnG00 -> 0
      EqnG11 -> 0

-- j == 0, v == -1
eqnBndV1 :: (Eq a, Floating a) => Eqn -> Coord a -> State a -> State a -> State a -> State a -> State a -> State a -> a
eqnBndV1 eqn c@(Coord u v) s@(State g q w) su _ _ suu _ =
  assert (v == -1) $
    case eqn of
      EqnG01 -> eqnG00 c s su suu
      EqnG22 -> q - 0
      EqnSW2 -> w - bndw c
      EqnG00 -> 0
      EqnG11 -> 0

eqnEndU1 :: Floating a => Eqn -> Coord a -> State a -> State a -> State a -> State a -> State a -> State a -> a
eqnEndU1 eqn c@(Coord u v) s@(State g q w) _ _ _ _ _ =
  case eqn of
    EqnG01 -> g - 1
    EqnG22 -> q - 0
    EqnSW2 -> w - bndw c
    EqnG00 -> 0
    EqnG11 -> 0

eqnEndV1 :: Floating a => Eqn -> Coord a -> State a -> State a -> State a -> State a -> State a -> State a -> a
eqnEndV1 eqn c@(Coord u v) (State g q w) _ _ _ _ _ =
  case eqn of
    EqnG01 -> g - 1
    EqnG22 -> q - 0
    EqnSW2 -> w - bndw c
    EqnG00 -> 0
    EqnG11 -> 0

equation1 :: (Eq a, Floating a) => Coord Int -> Eqn -> Coord a -> State a -> State a -> State a -> State a -> State a -> State a -> a
equation1 (Coord i j) =
  if
    | i == 0 && j == 0 -> eqnOrig1
    | i == 0 && j == np - 1 -> eqnEndU1
    | i == np - 1 && j == 0 -> eqnEndV1
    | i == 0 -> eqnBndU1
    | j == 0 -> eqnBndV1
    | i == j -> eqnAxis1
    | True -> eqnInt1

--------------------------------------------------------------------------------

g2s :: Storable a => Grid a -> State (Grid a)
g2s (Grid xs) =
  let [g, q, w] = L.takesV (replicate 3 (np * np)) xs
   in State (Grid g) (Grid q) (Grid w)

s2g :: Storable a => State (Grid a) -> Grid a
s2g (State (Grid g) (Grid q) (Grid w)) = Grid (L.vjoin [g, q, w])

eqn1 ::
  forall a.
  (Eq a, L.Field a, U.Unbox a) =>
  Eqn ->
  Coord (Grid a) ->
  State (Grid a) ->
  Grid a
eqn1 eqn c s =
  Grid $
    L.fromList
      [ eqnRow i j
        | i <- [0 .. np - 1],
          j <- [0 .. np - 1]
      ]
  where
    su = (derivU #>) <$> s
    sv = (derivV #>) <$> s
    suv = (derivUV #>) <$> s
    suu = (derivUU #>) <$> s
    svv = (derivVV #>) <$> s
    eqnRow :: Int -> Int -> a
    eqnRow i j =
      let at :: Grid a -> a
          at x = x ! (i, j)
       in equation1
            (Coord i j)
            eqn
            (at <$> c)
            (at <$> s)
            (at <$> su)
            (at <$> sv)
            (at <$> suv)
            (at <$> suu)
            (at <$> svv)

eqns ::
  forall a.
  (Eq a, L.Field a, U.Unbox a) =>
  Coord (Grid a) ->
  State (Grid a) ->
  Grid a
eqns c s =
  Grid $
    L.fromList
      [ eqnRow (i, j, eqn)
        | eqn <- [EqnG01, EqnG22, EqnSW2],
          i <- [0 .. np - 1],
          j <- [0 .. np - 1]
      ]
  where
    su = (derivU #>) <$> s
    sv = (derivV #>) <$> s
    suv = (derivUV #>) <$> s
    suu = (derivUU #>) <$> s
    svv = (derivVV #>) <$> s
    eqnRow :: (Int, Int, Eqn) -> a
    eqnRow (i, j, eqn) =
      let at :: Grid a -> a
          at x = x ! (i, j)
       in equation1
            (Coord i j)
            eqn
            (at <$> c)
            (at <$> s)
            (at <$> su)
            (at <$> sv)
            (at <$> suv)
            (at <$> suu)
            (at <$> svv)

eqns5 ::
  forall a.
  (Eq a, L.Field a, U.Unbox a) =>
  Coord (Grid a) ->
  State (Grid a) ->
  Grid a
eqns5 c s =
  Grid $
    L.fromList
      [ eqnRow (i, j, eqn)
        | eqn <- [EqnG01, EqnG22, EqnSW2, EqnG00, EqnG11],
          i <- [0 .. np - 1],
          j <- [0 .. np - 1]
      ]
  where
    su = (derivU #>) <$> s
    sv = (derivV #>) <$> s
    suv = (derivUV #>) <$> s
    suu = (derivUU #>) <$> s
    svv = (derivVV #>) <$> s
    eqnRow :: (Int, Int, Eqn) -> a
    eqnRow (i, j, eqn) =
      let at :: Grid a -> a
          at x = x ! (i, j)
       in equation1
            (Coord i j)
            eqn
            (at <$> c)
            (at <$> s)
            (at <$> su)
            (at <$> sv)
            (at <$> suv)
            (at <$> suu)
            (at <$> svv)

op ::
  forall a.
  (L.Container L.Vector a, Eq a, L.Field a, U.Unbox a) =>
  Coord (Grid a) ->
  State (Grid a) ->
  Op a
op (Coord u v) s =
  Op $
    L.fromRows
      [ jacRow (i, j, eqn)
        | eqn <- [EqnG01, EqnG22, EqnSW2],
          i <- [0 .. np - 1],
          j <- [0 .. np - 1]
      ]
  where
    su = (derivU #>) <$> s
    sv = (derivV #>) <$> s
    suv = (derivUV #>) <$> s
    suu = (derivUU #>) <$> s
    svv = (derivVV #>) <$> s
    jacRow :: (Int, Int, Eqn) -> L.Vector a
    jacRow (i, j, eqn) =
      L.vjoin $
        [ let con :: Grid a -> Dual a
              der :: Grid a -> Dual a
              con x = Dual (x ! (i, j)) 0
              der x = Dual (x ! (i, j)) 1
              scon :: State (Grid a) -> State (Dual a)
              sder :: State (Grid a) -> State (Dual a)
              scon (State g q w) = State (con g) (con q) (con w)
              sder (State g q w) = case var of
                VarG -> State (der g) (con q) (con w)
                VarQ -> State (con g) (der q) (con w)
                VarW -> State (con g) (con q) (der w)
           in ( ( dual $
                    equation1
                      (Coord i j)
                      eqn
                      (Coord (con u) (con v))
                      (sder s)
                      (scon su)
                      (scon sv)
                      (scon suv)
                      (scon suu)
                      (scon svv)
                )
                  `L.scale` (getOp delta L.! (i * np + j))
              )
                `L.add` ( ( dual $
                              equation1
                                (Coord i j)
                                eqn
                                (Coord (con u) (con v))
                                (scon s)
                                (sder su)
                                (scon sv)
                                (scon suv)
                                (scon suu)
                                (scon svv)
                          )
                            `L.scale` (getOp derivU L.! (i * np + j))
                        )
                `L.add` ( ( dual $
                              equation1
                                (Coord i j)
                                eqn
                                (Coord (con u) (con v))
                                (scon s)
                                (scon su)
                                (sder sv)
                                (scon suv)
                                (scon suu)
                                (scon svv)
                          )
                            `L.scale` (getOp derivV L.! (i * np + j))
                        )
                `L.add` ( ( dual $
                              equation1
                                (Coord i j)
                                eqn
                                (Coord (con u) (con v))
                                (scon s)
                                (scon su)
                                (scon sv)
                                (sder suv)
                                (scon suu)
                                (scon svv)
                          )
                            `L.scale` (getOp derivUV L.! (i * np + j))
                        )
                `L.add` ( ( dual $
                              equation1
                                (Coord i j)
                                eqn
                                (Coord (con u) (con v))
                                (scon s)
                                (scon su)
                                (scon sv)
                                (scon suv)
                                (sder suu)
                                (scon svv)
                          )
                            `L.scale` (getOp derivUU L.! (i * np + j))
                        )
                `L.add` ( ( dual $
                              equation1
                                (Coord i j)
                                eqn
                                (Coord (con u) (con v))
                                (scon s)
                                (scon su)
                                (scon sv)
                                (scon suv)
                                (scon suu)
                                (sder svv)
                          )
                            `L.scale` (getOp derivVV L.! (i * np + j))
                        )
          | var <- [VarG, VarQ, VarW]
        ]

op5 ::
  forall a.
  (L.Container L.Vector a, Eq a, L.Field a, U.Unbox a) =>
  Coord (Grid a) ->
  State (Grid a) ->
  Op a
op5 (Coord u v) s =
  Op $
    L.fromRows
      [ jacRow (i, j, eqn)
        | eqn <- [EqnG01, EqnG22, EqnSW2, EqnG00, EqnG11],
          i <- [0 .. np - 1],
          j <- [0 .. np - 1]
      ]
  where
    su = (derivU #>) <$> s
    sv = (derivV #>) <$> s
    suv = (derivUV #>) <$> s
    suu = (derivUU #>) <$> s
    svv = (derivVV #>) <$> s
    jacRow :: (Int, Int, Eqn) -> L.Vector a
    jacRow (i, j, eqn) =
      L.vjoin $
        [ let con :: Grid a -> Dual a
              der :: Grid a -> Dual a
              con x = Dual (x ! (i, j)) 0
              der x = Dual (x ! (i, j)) 1
              scon :: State (Grid a) -> State (Dual a)
              sder :: State (Grid a) -> State (Dual a)
              scon (State g q w) = State (con g) (con q) (con w)
              sder (State g q w) = case var of
                VarG -> State (der g) (con q) (con w)
                VarQ -> State (con g) (der q) (con w)
                VarW -> State (con g) (con q) (der w)
           in ( ( dual $
                    equation1
                      (Coord i j)
                      eqn
                      (Coord (con u) (con v))
                      (sder s)
                      (scon su)
                      (scon sv)
                      (scon suv)
                      (scon suu)
                      (scon svv)
                )
                  `L.scale` (getOp delta L.! (i * np + j))
              )
                `L.add` ( ( dual $
                              equation1
                                (Coord i j)
                                eqn
                                (Coord (con u) (con v))
                                (scon s)
                                (sder su)
                                (scon sv)
                                (scon suv)
                                (scon suu)
                                (scon svv)
                          )
                            `L.scale` (getOp derivU L.! (i * np + j))
                        )
                `L.add` ( ( dual $
                              equation1
                                (Coord i j)
                                eqn
                                (Coord (con u) (con v))
                                (scon s)
                                (scon su)
                                (sder sv)
                                (scon suv)
                                (scon suu)
                                (scon svv)
                          )
                            `L.scale` (getOp derivV L.! (i * np + j))
                        )
                `L.add` ( ( dual $
                              equation1
                                (Coord i j)
                                eqn
                                (Coord (con u) (con v))
                                (scon s)
                                (scon su)
                                (scon sv)
                                (sder suv)
                                (scon suu)
                                (scon svv)
                          )
                            `L.scale` (getOp derivUV L.! (i * np + j))
                        )
                `L.add` ( ( dual $
                              equation1
                                (Coord i j)
                                eqn
                                (Coord (con u) (con v))
                                (scon s)
                                (scon su)
                                (scon sv)
                                (scon suv)
                                (sder suu)
                                (scon svv)
                          )
                            `L.scale` (getOp derivUU L.! (i * np + j))
                        )
                `L.add` ( ( dual $
                              equation1
                                (Coord i j)
                                eqn
                                (Coord (con u) (con v))
                                (scon s)
                                (scon su)
                                (scon sv)
                                (scon suv)
                                (scon suu)
                                (sder svv)
                          )
                            `L.scale` (getOp derivVV L.! (i * np + j))
                        )
          | var <- [VarG, VarQ, VarW]
        ]

--------------------------------------------------------------------------------

getTime :: IO Double
getTime = do
  t <- getCPUTime
  return $ fromInteger t / 1.0e+12

main :: IO ()
main =
  do
    putStrLn "Spherically Symmetric Scalar Field in double-Null coordinates"
    let coords = Coord coordU (coordV @Double)
    let ini = State (gmap g guess) (gmap q guess) (gmap w guess)
    vars0 <- return (s2g ini)
    let done (n, vars) =
          let res = eqns coords (g2s vars)
           in -- let res = eqns5 coords (g2s vars)
              gmaxabs res < 1.0e-12 || n >= 5
    let iter (n, vars) = do
          let res = eqns coords (g2s vars)
          -- let res = eqns5 coords (g2s vars)
          t <- getTime
          putStrLn $
            "iter " ++ show n ++ " time " ++ show (round t) ++ " sec |res| "
              ++ show (gmaxabs res)
          let jac = op coords (g2s vars)
          let dvars = solve jac res
          -- let jac = op5 coords (g2s vars)
          -- let dvars = solveLS jac res
          return (n + 1, vars ^-^ dvars)
    (n1, vars1) <- iterateUntilM done iter (0, vars0)
    putStrLn $ "var.g " ++ show (gmap approx $ g $ g2s vars1)
    putStrLn $ "var.q " ++ show (gmap approx $ q $ g2s vars1)
    putStrLn $ "var.w " ++ show (gmap approx $ w $ g2s vars1)
    putStrLn $ "res.G01 " ++ show (gmap approx $ eqn1 EqnG01 coords (g2s vars1))
    putStrLn $ "res.G22 " ++ show (gmap approx $ eqn1 EqnG22 coords (g2s vars1))
    putStrLn $ "res.SW2 " ++ show (gmap approx $ eqn1 EqnSW2 coords (g2s vars1))
    putStrLn $ "res.G00 " ++ show (gmap approx $ eqn1 EqnG00 coords (g2s vars1))
    putStrLn $ "res.G11 " ++ show (gmap approx $ eqn1 EqnG11 coords (g2s vars1))
    let res1 = eqns coords (g2s vars1)
    t <- getTime
    putStrLn $
      "iter " ++ show n1 ++ " time " ++ show (round t) ++ " sec |res| "
        ++ show (gmaxabs res1)
    let res5 = eqns5 coords (g2s vars1)
    putStrLn $ "|res5| " ++ show (gmaxabs res5)
    -- let h = 2 / (fromIntegral np - 1) :: Double
    -- let alpha = 0.01 * h ^ 3
    -- let smooth = delta - alpha *^ (derivUU * derivUU + derivVV * derivVV)
    -- let vars1s = s2g $ (expfilter #>) <$> g2s vars1
    -- let res5s = eqns5 coords (g2s vars1s)
    -- putStrLn $ "|res5s| " ++ show (gmaxabs res5s)
    putStrLn "Done."
  where
    approx :: RealFrac a => a -> a
    approx x = fromInteger (round (1.0e+10 * x)) / 1.0e+10
