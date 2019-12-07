{-# LANGUAGE TypeApplications #-}
{-# OPTIONS_GHC -fno-warn-type-defaults #-}

module DerivProps where

import Data.VectorSpace
import qualified GHC.Exts as GHC
import Numeric.LinearAlgebra hiding ((#>), (===), inv)
import SSSFN
import Test.QuickCheck

default (Int)

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

prop_derivU_coordU :: Property
prop_derivU_coordU = derivU #> coordU @Double ~~~ gpure 1

prop_derivV_coordV :: Property
prop_derivV_coordV = derivV #> coordV @Double ~~~ gpure 1

prop_derivU_coordV :: Property
prop_derivU_coordV = derivU #> coordV @Double ~~~ gpure 0

prop_derivV_coordU :: Property
prop_derivV_coordU = derivV #> coordU @Double ~~~ gpure 0

prop_derivT_coordT :: Property
prop_derivT_coordT = derivT #> coordT @Double ~~~ gpure 1

prop_derivR_coordR :: Property
prop_derivR_coordR = derivR #> coordR @Double ~~~ gpure 1

prop_derivT_coordR :: Property
prop_derivT_coordR = derivT #> coordR @Double ~~~ gpure 0

prop_derivR_coordT :: Property
prop_derivR_coordT = derivR #> coordT @Double ~~~ gpure 0

prop_derivU_coordR :: Property
prop_derivU_coordR = derivU #> coordR @Double ~~~ gpure (-1)

prop_derivV_coordR :: Property
prop_derivV_coordR = derivV #> coordR @Double ~~~ gpure 1

prop_symmetrizeR_coordT :: Property
prop_symmetrizeR_coordT = symmetrizeR #> coordT @Double ~~~ coordT

prop_symmetrizeR_coordR :: Property
prop_symmetrizeR_coordR = symmetrizeR #> coordR @Double ~~~ gpure 0

prop_zproject_const :: Property
prop_zproject_const = zproject #> gpure 1 ~~~ gpure (0 :: Double)

prop_zproject_coordR :: Property
prop_zproject_coordR = zproject #> (coordR :: Grid Double) ~~~ gpure 0

prop_zproject_coordR2 :: Property
prop_zproject_coordR2 = zproject #> (coordR2 :: Grid Double) ~~~ coordR2
  where
    coordR2 = gmap (^ 2) coordR

prop_zproject_coordR2T :: Property
prop_zproject_coordR2T = zproject #> (coordR2T :: Grid Double) ~~~ coordR2T
  where
    coordR2T = gzipWith (\t r -> t * r ^ 2) coordT coordR

prop_zproject_coordR4 :: Property
prop_zproject_coordR4 = zproject #> (coordR4 :: Grid Double) ~~~ coordR4
  where
    coordR4 = gmap (^ 4) coordR

prop_r2project_const :: Property
prop_r2project_const = r2project #> (gpure 1 :: Grid Double) ~~~ gpure 0

prop_r2project_coordR :: Property
prop_r2project_coordR = r2project #> (coordR :: Grid Double) ~~~ gpure 0

prop_r2project_coordR2 :: Property
prop_r2project_coordR2 = r2project #> (coordR2 :: Grid Double) ~~~ gpure 1
  where
    coordR2 = gmap (^ 2) coordR

prop_r2project_coordR2T :: Property
prop_r2project_coordR2T = r2project #> (coordR2T :: Grid Double) ~~~ coordT
  where
    coordR2T = gzipWith (\t r -> t * r ^ 2) coordT coordR

prop_r2project_coordR4 :: Property
prop_r2project_coordR4 = r2project #> (coordR4 :: Grid Double) ~~~ coordR2
  where
    coordR2 = gmap (^ 2) coordR
    coordR4 = gmap (^ 4) coordR

prop_r1derivR_coordR :: Property
prop_r1derivR_coordR = r1derivR #> (coordR :: Grid Double) ~~~ gpure 0

prop_r1derivR_coordR2 :: Property
prop_r1derivR_coordR2 = r1derivR #> (coordR2 :: Grid Double) ~~~ gpure 2
  where
    coordR2 = gmap (^ 2) coordR

prop_r1derivR_coordR2T :: Property
prop_r1derivR_coordR2T = r1derivR #> (coordR2T :: Grid Double) ~~~ 2 *^ coordT
  where
    coordR2T = gzipWith (\t r -> t * r ^ 2) coordT coordR

prop_r1derivR_coordR4 :: Property
prop_r1derivR_coordR4 = r1derivR #> (coordR4 :: Grid Double) ~~~ 4 *^ coordR2
  where
    coordR2 = gmap (^ 2) coordR
    coordR4 = gmap (^ 4) coordR

prop_r1derivU_coordR :: Property
prop_r1derivU_coordR = r1derivU #> (coordR :: Grid Double) ~~~ gpure 0

prop_r1derivU_coordR2 :: Property
prop_r1derivU_coordR2 = r1derivU #> (coordR2 :: Grid Double) ~~~ gpure (-2)
  where
    coordR2 = gmap (^ 2) coordR

prop_r1derivU_coordR4 :: Property
prop_r1derivU_coordR4 = r1derivU #> (coordR4 :: Grid Double) ~~~ (-4) *^ coordR2
  where
    coordR2 = gmap (^ 2) coordR
    coordR4 = gmap (^ 4) coordR

prop_r1derivV_coordR :: Property
prop_r1derivV_coordR = r1derivV #> (coordR :: Grid Double) ~~~ gpure 0

prop_r1derivV_coordR2 :: Property
prop_r1derivV_coordR2 = r1derivV #> (coordR2 :: Grid Double) ~~~ gpure 2
  where
    coordR2 = gmap (^ 2) coordR

prop_r1derivV_coordR4 :: Property
prop_r1derivV_coordR4 = r1derivV #> (coordR4 :: Grid Double) ~~~ 4 *^ coordR2
  where
    coordR2 = gmap (^ 2) coordR
    coordR4 = gmap (^ 4) coordR
