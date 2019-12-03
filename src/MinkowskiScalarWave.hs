{-# LANGUAGE TypeApplications #-}

module MinkowskiScalarWave where

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

potential :: forall a. (L.Container L.Vector a, Fractional a) => Grid a
potential = gzipWith (\u v -> 0) (coordU @a) (coordV @a)

incoming :: forall a. (L.Container L.Vector a, Floating a) => Grid a
incoming = gsample bcond

--------------------------------------------------------------------------------

op :: (L.Field a, Ord a, U.Unbox a) => Op a
op = intr * dalembert + bndu * incomingU + bndv * incomingV + orig * originUV
  where
    bndu = incomingU - originUV
    bndv = incomingV - originUV
    orig = originUV
    intr = 1 - (bndu + bndv + orig)

rhs :: (Floating a, L.Numeric a, Ord a) => Grid a
rhs =
  intr #> potential
    ^+^ bndu #> incoming
    ^+^ bndv #> incoming
    ^+^ orig #> incoming
  where
    bndu = incomingU - originUV
    bndv = incomingV - originUV
    orig = originUV
    intr = 1 - (bndu + bndv + orig)

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
    putStrLn "1+1D Scalar Field in the Minkowski spacetime in double-Null coordinates"
    -- putStrLn $ "intr " ++ show (intr :: Op Double)
    -- putStrLn $ "bndu " ++ show (bndu :: Op Double)
    -- putStrLn $ "bndv " ++ show (bndv :: Op Double)
    -- putStrLn $ "orig " ++ show (orig :: Op Double)
    -- putStrLn $ "dal " ++ show (omap approx (dalembert @Double))
    -- putStrLn $ "op " ++ show (omap approx (op @Double))
    -- let inc = incoming @Double
    -- putStrLn $ "inc " ++ show (gmap approx inc)
    -- putStrLn $ "rhs " ++ show (gmap approx (rhs @Double))
    -- let ires = dalembert #> inc
    -- putStrLn $ "ires " ++ show (gmap approx ires)
    -- putStrLn $ "|ires| " ++ show (gmaxabs ires)
    -- let jres = op #> inc ^-^ rhs
    -- putStrLn $ "jres " ++ show (gmap approx jres)
    -- putStrLn $ "|jres| " ++ show (gmaxabs jres)

    let sol = solve op (rhs @Double)
    putStrLn $ "sol " ++ show (gmap approx sol)
    let res = op #> sol ^-^ rhs
    putStrLn $ "res " ++ show (gmap approx res)
    putStrLn $ "|res| " ++ show (gmaxabs res)
    let err = sol ^-^ incoming
    putStrLn $ "err " ++ show (gmap approx err)
    putStrLn $ "|err| " ++ show (gmaxabs err)
    putStrLn "Done."
  where
    approx :: RealFrac a => a -> a
    approx x = fromInteger (round (1.0e+4 * x)) / 1.0e+4
-- intr = 1 - (incomingU + incomingV - originUV)
-- bndu = incomingU - originUV
-- bndv = incomingV - originUV
-- orig = originUV

--------------------------------------------------------------------------------

-- quantum measurements: the quantum system measures the state of the
-- classical system

-- Beatrice Bonga: conserved current from varying boundary term from
-- lagrangian <https://arxiv.org/abs/1911.04514>
