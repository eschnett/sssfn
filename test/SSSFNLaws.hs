{-# LANGUAGE TypeApplications #-}

module SSSFNLaws where

import Data.VectorSpace
import qualified GHC.Exts as GHC
import Numeric.LinearAlgebra hiding ((#>), (===), inv)
import SSSFN
import Test.QuickCheck

--------------------------------------------------------------------------------

infix 4 ~~~

(~~~) ::
  (GHC.IsList as, Show as, a ~ GHC.Item as, RealFloat a, Show a) =>
  as ->
  as ->
  Property
xs ~~~ ys =
  let xs' = GHC.toList xs
      ys' = GHC.toList ys
      sc = 1 `max` maximum (map abs xs') `max` maximum (map abs ys')
      di = maximum (zipWith absdiff xs' ys')
   in counterexample
        ( show xs
            ++ " ~~~ "
            ++ show ys
            ++ " (difference "
            ++ show di
            ++ ", scale "
            ++ show sc
            ++ ")"
        )
        (length xs' == length ys' && di <= eps * sc)
  where
    absdiff x y = abs (x - y)
    eps = 100 * peps

--------------------------------------------------------------------------------

prop_Grid_add_assoc :: Grid Double -> Grid Double -> Grid Double -> Property
prop_Grid_add_assoc x y z = (x ^+^ y) ^+^ z ~~~ x ^+^ (y ^+^ z)

prop_Grid_add_lunit :: Grid Double -> Property
prop_Grid_add_lunit x = zeroV ^+^ x ~~~ x

prop_Grid_add_runit :: Grid Double -> Property
prop_Grid_add_runit x = x ~~~ x ^+^ zeroV

prop_Grid_add_comm :: Grid Double -> Grid Double -> Property
prop_Grid_add_comm x y = x ^+^ y ~~~ y ^+^ x

prop_Grid_add_linv :: Grid Double -> Property
prop_Grid_add_linv x = negateV x ^+^ x ~~~ zeroV

prop_Grid_add_rinv :: Grid Double -> Property
prop_Grid_add_rinv x = x ^+^ negateV x ~~~ zeroV

prop_Grid_add_scale_dist :: Double -> Grid Double -> Grid Double -> Property
prop_Grid_add_scale_dist a x y = a *^ (x ^+^ y) ~~~ a *^ x ^+^ a *^ y

prop_Grid_scale_assoc :: Double -> Double -> Grid Double -> Property
prop_Grid_scale_assoc a b x = a *^ (b *^ x) ~~~ (a * b) *^ x

prop_Grid_scale_lunit :: Grid Double -> Property
prop_Grid_scale_lunit x = 1 *^ x ~~~ x

prop_Grid_scale_szero :: Grid Double -> Property
prop_Grid_scale_szero x = 0 *^ x ~~~ zeroV

prop_Grid_scale_zero :: Double -> Property
prop_Grid_scale_zero a = a *^ (zeroV :: Grid Double) ~~~ zeroV

prop_Grid_scale_comm :: Double -> Double -> Grid Double -> Property
prop_Grid_scale_comm a b x = a *^ (b *^ x) ~~~ b *^ (a *^ x)

prop_Grid_scale_Grid_add_dist :: Double -> Double -> Grid Double -> Property
prop_Grid_scale_Grid_add_dist a b x = (a + b) *^ x ~~~ a *^ x ^+^ b *^ x

--------------------------------------------------------------------------------

prop_Op_add_assoc :: Op Double -> Op Double -> Op Double -> Property
prop_Op_add_assoc x y z = (x + y) + z ~~~ x + (y + z)

prop_Op_add_lunit :: Op Double -> Property
prop_Op_add_lunit x = 0 + x ~~~ x

prop_Op_add_runit :: Op Double -> Property
prop_Op_add_runit x = x ~~~ x + 0

prop_Op_add_comm :: Op Double -> Op Double -> Property
prop_Op_add_comm x y = x + y ~~~ y + x

prop_Op_add_linv :: Op Double -> Property
prop_Op_add_linv x = negate x + x ~~~ 0

prop_Op_add_rinv :: Op Double -> Property
prop_Op_add_rinv x = x + negate x ~~~ 0

prop_Op_add_scale_dist :: Double -> Op Double -> Op Double -> Property
prop_Op_add_scale_dist a x y = a *^ (x + y) ~~~ a *^ x + a *^ y

prop_Op_mul_assoc :: Op Double -> Op Double -> Op Double -> Property
prop_Op_mul_assoc x y z = (x * y) * z ~~~ x * (y * z)

prop_Op_mul_lunit :: Op Double -> Property
prop_Op_mul_lunit x = 1 * x ~~~ x

prop_Op_mul_runit :: Op Double -> Property
prop_Op_mul_runit x = x ~~~ x * 1

prop_Op_mul_scale_assoc :: Double -> Op Double -> Op Double -> Property
prop_Op_mul_scale_assoc a x y = a *^ (x * y) ~~~ (a *^ x) * y

prop_Op_scale_assoc :: Double -> Double -> Op Double -> Property
prop_Op_scale_assoc a b x = a *^ (b *^ x) ~~~ (a * b) *^ x

prop_Op_scale_lunit :: Op Double -> Property
prop_Op_scale_lunit x = 1 *^ x ~~~ x

prop_Op_scale_lzero :: Op Double -> Property
prop_Op_scale_lzero x = 0 *^ x ~~~ 0

prop_Op_scale_rzero :: Double -> Property
prop_Op_scale_rzero a = a *^ (0 :: Op Double) ~~~ 0

prop_Op_scale_comm :: Double -> Double -> Op Double -> Property
prop_Op_scale_comm a b x = a *^ (b *^ x) ~~~ b *^ (a *^ x)

prop_Op_scale_add_dist :: Double -> Double -> Op Double -> Property
prop_Op_scale_add_dist a b x = (a + b) *^ x ~~~ (a *^ x) + (b *^ x)

--------------------------------------------------------------------------------

-- prop_deriv_coord :: Property
-- prop_deriv_coord =
--   derivs L.#> coords @Double ~~~ cmap (const 1) (coords @Double)

prop_derivU_coordU :: Property
prop_derivU_coordU = derivU #> coordU @Double ~~~ gpure 1

prop_derivU_coordV :: Property
prop_derivU_coordV = derivU #> coordV @Double ~~~ gpure 0

prop_derivV_coordU :: Property
prop_derivV_coordU = derivV #> coordU @Double ~~~ gpure 0

prop_derivV_coordV :: Property
prop_derivV_coordV = derivV #> coordV @Double ~~~ gpure 1
-- -- W D + (W D)^T = W B
-- prop_derivU_bndU :: Property
-- prop_derivU_bndU =
--   metric * derivU + transpose (metric * derivU) ~~~ metric * bndU @Double
