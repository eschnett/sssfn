{-# LANGUAGE TypeApplications #-}
{-# OPTIONS_GHC -fno-warn-type-defaults #-}

module ScalarWave where

import Control.Exception (assert)
import qualified Data.Vector.Unboxed as U
import Data.VectorSpace
import Dual
import Equations
import Foreign
import qualified Numeric.LinearAlgebra as L
import SSSFN

default (Int)

--------------------------------------------------------------------------------

bndw :: Floating a => Coord a -> a
bndw (Coord u v) = swamp * cos (pi / 4 * t) * cos (pi / 4 * r)
  where
    t = v + u
    r = v - u
    swamp = 0.01

guess1 :: Floating a => Coord a -> State a
-- guess1 c = State 1 0 (bndw c)
guess1 c = State 1 0 0

guess :: (L.Container L.Vector a, Floating a) => Grid (State a)
guess = gsample (\(u, v) -> guess1 (Coord u v))

--------------------------------------------------------------------------------

eqn1 :: Fractional a => Eqn -> Coord a -> State a -> State a -> State a -> State a -> State a -> State a -> a
eqn1 eqn c s su sv suv suu svv =
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
    | True -> eqn1

--------------------------------------------------------------------------------

g2s :: Storable a => Grid a -> State (Grid a)
g2s (Grid xs) =
  let [g, q, w] = L.takesV (replicate 3 (np * np)) xs
   in State (Grid g) (Grid q) (Grid w)

s2g :: Storable a => State (Grid a) -> Grid a
s2g (State (Grid g) (Grid q) (Grid w)) = Grid (L.vjoin [g, q, w])

eqns ::
  forall a.
  (Eq a, L.Field a, U.Unbox a) =>
  Coord (Grid a) ->
  State (Grid a) ->
  Grid a
eqns c s =
  Grid $
    L.fromList
      [ resRow (i, j, eqn)
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
    resRow :: (Int, Int, Eqn) -> a
    resRow (i, j, eqn) =
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

--------------------------------------------------------------------------------

main :: IO ()
main =
  do
    putStrLn "Spherically Symmetric Scalar Field in double-Null coordinates"
    let coords = Coord coordU (coordV @Double)
    let ini = State (gmap g guess) (gmap q guess) (gmap w guess)
    vars0 <- return (s2g ini)
    let iter vars = do
          -- putStrLn $ "var.g " ++ show (gmap approx $ g $ g2s vars)
          -- putStrLn $ "var.q " ++ show (gmap approx $ q $ g2s vars)
          -- putStrLn $ "var.w " ++ show (gmap approx $ w $ g2s vars)
          let res = eqns coords (g2s vars)
          -- putStrLn $ "res.G01 " ++ show (gmap approx $ g $ g2s res)
          -- putStrLn $ "res.G22 " ++ show (gmap approx $ q $ g2s res)
          -- putStrLn $ "res.SW2 " ++ show (gmap approx $ w $ g2s res)
          putStrLn $ "|res| " ++ show (gmaxabs res)
          let jac = op coords (g2s vars)
          let dvars = solve jac res
          return (vars ^-^ dvars)
    vars1 <- iter vars0
    vars2 <- iter vars1
    vars3 <- iter vars2
    vars4 <- iter vars3
    putStrLn "Done."
  where
    approx :: RealFrac a => a -> a
    approx x = fromInteger (round (1.0e+4 * x)) / 1.0e+4
