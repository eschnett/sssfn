{-# LANGUAGE TypeApplications #-}

import qualified Data.Vector.Unboxed as U
import Data.VectorSpace
import qualified Numeric.LinearAlgebra as L
import SSSFN

--------------------------------------------------------------------------------

-- | Scalar wave in double-null coordinates:
-- > (- dt^2 + dx^2) u = 0
-- > t = v + u   dt = dv + du
-- > x = v - u   dx = dv - du
-- > (- (dv + du)^2 + (dv - du)^2) u = 0
-- > - 4 dv du u = 0
dalembert :: (Eq a, L.Field a, U.Unbox a) => Op a
dalembert = -4 *^ (derivV * derivU)

incomingU :: forall a. (L.Container L.Vector a, Fractional a, Ord a) => Op a
incomingU = omap isneg (bndU @a)
  where
    isneg a = fromIntegral (fromEnum (a < 0))

incomingV :: forall a. (L.Container L.Vector a, Fractional a, Ord a) => Op a
incomingV = omap isneg (bndV @a)
  where
    isneg a = fromIntegral (fromEnum (a < 0))

originUV :: forall a. (L.Container L.Vector a, Fractional a, Ord a) => Op a
originUV = ozipWith (*) incomingU incomingV

potential ::
  forall a.
  (L.Container L.Vector a, Fractional a, L.Indexable (L.Vector a) a) =>
  Grid a
potential = gzipWith (\u v -> 0) (coordU @a) (coordV @a)

incoming ::
  forall a.
  (L.Container L.Vector a, Floating a, L.Indexable (L.Vector a) a) =>
  Grid a
incoming = gsample bcond

--------------------------------------------------------------------------------

op :: (L.Field a, Ord a, U.Unbox a) => Op a
op =
  intr * dalembert * intr
    + bndu * incomingU * bndu
    + bndv * incomingV * bndv
    + orig * originUV * orig
  where
    intr = 1 - (incomingU + incomingV - originUV)
    bndu = incomingU - originUV
    bndv = incomingV - originUV
    orig = originUV

rhs :: (Floating a, L.Indexable (L.Vector a) a, L.Numeric a, Ord a) => Grid a
rhs =
  intr #> potential
    ^+^ bndu #> incoming
    ^+^ bndv #> incoming
    ^+^ orig #> incoming
  where
    intr = 1 - (incomingU + incomingV - originUV)
    bndu = incomingU - originUV
    bndv = incomingV - originUV
    orig = originUV

--------------------------------------------------------------------------------

bcond :: Floating a => (a, a) -> a
bcond (u, v) = cos (pi / 4 * t) * cos (pi / 4 * x)
  where
    t = v + u
    x = v - u

--------------------------------------------------------------------------------

main :: IO ()
main =
  do
    putStrLn "Spherically Symmetric Scalar Field in double-Null coordinates"
    putStrLn $ "u " ++ show (coordU @Double)
    putStrLn $ "v " ++ show (coordV @Double)
    putStrLn $ "inc " ++ show (gmap approx (incoming @Double))
    let res = dalembert #> incoming @Double
    putStrLn $ "res " ++ show (gmap approx res)
    putStrLn $ "|res| " ++ show (gmaxabs res)
    -- -- putStrLn $ "dal " ++ show (dalembert @Double)
    -- -- putStrLn $ "incu " ++ show (incomingU @Double)
    -- -- putStrLn $ "incv " ++ show (incomingV @Double)
    -- -- putStrLn $ "orig " ++ show (originUV @Double)
    -- -- putStrLn $ "op " ++ show (op @Double)
    -- putStrLn $ "pot " ++ show (potential @Double)
    -- putStrLn $ "rhs " ++ show (rhs @Double)
    -- let sol = solve op (rhs @Double)
    -- putStrLn $ "sol " ++ show sol
    -- let res = dalembert #> sol
    -- putStrLn $ "res " ++ show res
    -- putStrLn $ "|res| " ++ show (gmaxabs res)
    putStrLn "Done."
  where
    approx :: RealFrac a => a -> a
    approx x = fromInteger (round (1.0e+4 * x)) / 1.0e+4
