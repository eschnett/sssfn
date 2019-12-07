{-# LANGUAGE TypeApplications #-}
{-# OPTIONS_GHC -fno-warn-type-defaults #-}

module SSSFN
  ( -- * Domain
    np,

    -- * Grid
    Grid (..),
    gunit,
    (!),
    gcoerce,

    -- * Operators
    Op (..),
    ounit,
    ocoerce,

    -- * Grids are functors etc.
    gmap,
    gpure,
    gzipWith,
    gfoldMap,
    gcount,
    gmaxabs,
    gsum,

    -- * Operators are functors etc.
    omap,
    opure,
    ozipWith,
    ofoldMap,
    ocount,
    omaxabs,
    osum,

    -- * Operators
    (#>),
    transpose,
    inv,

    -- * Coordinates
    isOrigin,
    coordU,
    coordV,
    coordT,
    coordR,
    coordR1, -- remove?

    -- * Sampling functions
    gsample,

    -- * Evaluating functions
    geval,
    gevalCoeffs,

    -- * Derivatives
    vandermonde,
    expfilter,
    delta,
    symmetrizeU,
    symmetrizeV,
    symmetrizeR,
    derivU,
    derivV,
    derivUU,
    derivUV,
    derivVV,
    derivT,
    derivR,
    derivTT,
    derivRR,
    r1derivR,
    r1derivU,
    r1derivV,
    r2project,
    zproject,
    bndU,
    bndV,
    metric,

    -- * Systems of equations
    solve,
    solveLS,

    -- * MaxAbs
    maxabs,
    getMaxAbs,
  )
where

import Control.Exception (assert)
import Data.Maybe
import Data.Poly hiding (scale)
import Data.Semigroup hiding ((<>))
import qualified Data.Vector.Unboxed as U
import Data.VectorSpace
import Foreign
import qualified GHC.Exts as GHC
import Numeric.LinearAlgebra hiding ((!), (#>), inv)
import qualified Numeric.LinearAlgebra as L
import Numeric.LinearAlgebra.Devel
import qualified Test.QuickCheck as QC
import Prelude hiding ((<>))
import qualified Prelude as P

default (Int)

--------------------------------------------------------------------------------

-- | Grid size
np :: Int
np = 5

--------------------------------------------------------------------------------

-- | Grid datatype
newtype Grid a = Grid {getGrid :: Vector a}
  deriving (Eq, Ord, Read, Show)

instance (QC.Arbitrary a, Storable a) => QC.Arbitrary (Grid a) where

  arbitrary = Grid . fromList <$> QC.vector (np ^ 2)

  shrink = map (Grid . fromList) . traverse QC.shrink . toList . getGrid

instance (Storable a) => GHC.IsList (Grid a) where

  type Item (Grid a) = a

  fromList = Grid . fromList

  toList = toList . getGrid

(!) :: (Container Vector a) => Grid a -> (Int, Int) -> a
(!) (Grid x) (i, j) = x `atIndex` (i * np + j)

gunit :: (Num a, Storable a) => (Int, Int) -> Grid a
gunit (i, j) =
  let ij = i * np + j
   in Grid $ fromList [if k == ij then 1 else 0 | k <- [0 .. np ^ 2 -1]]

gcoerce :: forall a b. (Storable a, Storable b) => Grid a -> Grid b
gcoerce (Grid xss) =
  assert (na * sa `mod` sb == 0) $
    (Grid $ fromByteString $ toByteString xss)
  where
    na = cols (asRow xss)
    sa = sizeOf @a undefined
    sb = sizeOf @b undefined

--------------------------------------------------------------------------------

-- | Operator datatype
newtype Op a = Op {getOp :: Matrix a}
  deriving (Read, Show)

instance (QC.Arbitrary a, Element a) => QC.Arbitrary (Op a) where

  arbitrary = Op . reshape (np ^ 2) . fromList <$> QC.vector (np ^ 4)

  shrink =
    map (Op . reshape (np ^ 2) . fromList)
      . traverse QC.shrink
      . (toList . flatten . getOp)

instance Element a => GHC.IsList (Op a) where

  type Item (Op a) = a

  fromList = Op . reshape (np ^ 2) . fromList

  toList = toList . flatten . getOp

ounit :: (Element a, Num a) => Op a
ounit = Op $ ident (np ^ 2)

ocoerce :: forall a b. (Container Vector a, Num a, Storable b) => Op a -> Op b
ocoerce (Op x) =
  assert (ma * sa `mod` sb == 0 && na * sa `mod` sb == 0) $
    (Op $ reshape nb $ getGrid $ gcoerce $ Grid $ flatten x)
  where
    (ma, na) = size x
    sa = sizeOf @a undefined
    sb = sizeOf @b undefined
    nb = na * sa `div` sb

--------------------------------------------------------------------------------

-- | Grids are functors etc.
-- gmap :: (Container Vector a, Element b) => (a -> b) -> Grid a -> Grid b
-- gmap f (Grid x) = Grid (cmap f x)
gmap :: (Storable a, Storable b) => (a -> b) -> Grid a -> Grid b
gmap f (Grid x) = Grid (mapVectorWithIndex (\i a -> f a) x)

gpure :: Container Vector a => a -> Grid a
gpure a = Grid $ konst a (np ^ 2)

gzipWith ::
  (Storable a, Storable b, Storable c) =>
  (a -> b -> c) ->
  Grid a ->
  Grid b ->
  Grid c
gzipWith f (Grid x) (Grid y) = Grid (zipVectorWith f x y)

gfoldMap :: (Storable a, Monoid b) => (a -> b) -> Grid a -> b
gfoldMap f (Grid x) = foldVector (\a b -> f a P.<> b) mempty x

gcount :: Container Vector a => Grid a -> Int
gcount (Grid x) = size x

gmaxabs :: (Container Vector a, Num a, Ord a) => Grid a -> a
gmaxabs = getMaxAbs . gfoldMap maxabs

gsum :: (Container Vector a) => Grid a -> a
gsum (Grid x) = sumElements x

--------------------------------------------------------------------------------

-- | Operators are functors etc.
-- omap :: (Container Vector a, Num a, Element b) => (a -> b) -> Op a -> Op b
-- omap f (Op x) = Op (cmap f x)
omap :: (L.Element a, Storable b) => (a -> b) -> Op a -> Op b
omap f (Op x) = Op $ reshape n $ getGrid $ gmap f $ Grid $ flatten x
  where
    n = cols x

opure :: (Container Vector a, Num a) => a -> Op a
opure a = Op (konst a (np ^ 2, np ^ 2))

ozipWith ::
  (Element a, Element b, Storable c) =>
  (a -> b -> c) ->
  Op a ->
  Op b ->
  Op c
ozipWith f (Op x) (Op y) =
  Op $ reshape (np ^ 2) $ zipVectorWith f (flatten x) (flatten y)

ofoldMap :: (Element a, Monoid b) => (a -> b) -> Op a -> b
ofoldMap f (Op x) = foldVector (\a b -> f a P.<> b) mempty (flatten x)

ocount :: (Container Vector a, Num a) => Op a -> Int
ocount (Op x) = let (ni, nj) = size x in ni * nj

omaxabs :: (Container Vector a, Num a, Ord a) => Op a -> a
omaxabs = getMaxAbs . ofoldMap maxabs

osum :: (Container Vector a, Num a) => Op a -> a
osum (Op x) = sumElements x

--------------------------------------------------------------------------------

-- | Grids are an additive group
instance
  (Container Vector a, Num a) =>
  AdditiveGroup (Grid a)
  where

  zeroV = gpure 0

  (^+^) = gzipWith (+)

  negateV = gmap negate

-- | Grids are a vector space
instance
  (Container Vector a, Num a) =>
  VectorSpace (Grid a)
  where

  type Scalar (Grid a) = a

  (*^) a = gmap (a *)

--------------------------------------------------------------------------------

-- | Operators form an additive group
instance (Container Vector a, Num a) => AdditiveGroup (Op a) where

  zeroV = opure 0

  (^+^) = ozipWith (+)

  negateV = omap negate

-- | Operators form a vector space
instance (Container Vector a, Num a) => VectorSpace (Op a) where

  type Scalar (Op a) = a

  (*^) a = omap (a *)

-- | Operators form a ring
instance (Container Vector a, Num a, Numeric a) => Num (Op a) where

  (+) = (^+^)

  negate = negateV

  (-) = (^-^)

  Op x * Op y = Op (x <> y)

  abs = omap abs

  signum = omap signum

  fromInteger i = fromInteger i *^ ounit

-- | Matrix-vector multiplication
infixr 8 #>

(#>) :: Numeric a => Op a -> Grid a -> Grid a
Op x #> Grid y = Grid (x L.#> y)

-- | Matrix transpose
transpose :: Transposable (Matrix a) (Matrix a) => Op a -> Op a
transpose (Op x) = Op (tr x)

-- | Matrix inverse
inv :: Field a => Op a -> Op a
inv (Op x) = Op (L.inv x)

--------------------------------------------------------------------------------

-- | Coordinate of a grid point
coord :: Floating a => Int -> a
coord i =
  let imin = 0
      imax = np - 1
      umin = -1
      umax = 1
      u =
        fromIntegral (i - imax) / fromIntegral (imin - imax) * umin
          + fromIntegral (i - imin) / fromIntegral (imax - imin) * umax
   in sin (pi / 2 * u)

-- | Coordinates of all grid points
coords :: (Floating a, Storable a) => Vector a
coords = fromList [coord i | i <- [0 .. np -1]]

poly2vector :: (Num a, Storable a, U.Unbox a) => UPoly a -> Vector a
poly2vector p =
  let v = unPoly p
      n = U.length v
   in fromList [if i < n then v U.! i else 0 | i <- [0 .. np -1]]

-- | Evalute a polynomial at a grid point
evalcoeffs :: (Eq a, Num a, Storable a, U.Unbox a) => a -> Vector a
evalcoeffs u = fromList [eval (mon i) u | i <- [0 .. np -1]]
  where
    mon :: (Eq a, Num a, U.Unbox a) => Int -> UPoly a
    mon i = monomial (fromIntegral i) 1

-- | Evaluate a polynomial at all grid points
gvandermonde :: (Element a, Eq a, Floating a, U.Unbox a) => Matrix a
gvandermonde = fromRows [evalcoeffs u | u <- toList coords]

--------------------------------------------------------------------------------

-- Note: hmatrix is row-major by default
mkGrid :: (Num a, Storable a) => (Int -> a) -> (Int -> a) -> Grid a
mkGrid fx fy =
  Grid $ fromList [fx i * fy j | i <- [0 .. np -1], j <- [0 .. np -1]]

isOrigin :: Grid Bool
isOrigin =
  Grid $
    fromList
      [ i == j
        | i <- [0 .. np -1],
          j <- [0 .. np -1]
      ]

-- | Coordinates of all grid points
coordU :: (Container Vector a, Floating a) => Grid a
coordU = mkGrid (coords `atIndex`) (const 1)

-- | Coordinates of all grid points
coordV :: (Container Vector a, Floating a) => Grid a
coordV = mkGrid (const 1) (coords `atIndex`)

-- | Combined coordinates
coordT :: (Container Vector a, Floating a) => Grid a
coordT = coordV ^+^ coordU

coordR :: (Container Vector a, Floating a) => Grid a
coordR = coordV ^-^ coordU

coordR1 :: (Container Vector a, RealFloat a) => Grid a
coordR1 = gzipWith (\isOrig r -> if isOrig then 0 else recip r) isOrigin coordR

--------------------------------------------------------------------------------

-- | Sample a function at grid points
gsample ::
  (Container Vector a, Floating a, Storable b) =>
  ((a, a) -> b) ->
  Grid b
gsample f = gzipWith (curry f) (coordU) (coordV)

--------------------------------------------------------------------------------

-- | Evaluate a grid at a point
-- f / Identity
geval ::
  ( Container Vector a,
    Eq a,
    Field a,
    U.Unbox a,
    Container Vector b,
    VectorSpace b,
    Scalar b ~ a
  ) =>
  Grid b ->
  (a, a) ->
  b
geval g (u, v) =
  let ucoeffs = tr (L.inv gvandermonde) L.#> evalcoeffs u
      vcoeffs = tr (L.inv gvandermonde) L.#> evalcoeffs v
   in sumV
        [ ((ucoeffs `atIndex` i) * (vcoeffs `atIndex` j)) *^ (g ! (i, j))
          | i <- [0 .. np - 1],
            j <- [0 .. np - 1]
        ]

gevalCoeffs ::
  ( Container Vector a,
    Eq a,
    Field a,
    U.Unbox a
  ) =>
  (a, a) ->
  Grid a
gevalCoeffs (u, v) =
  let ucoeffs = tr (L.inv gvandermonde) L.#> evalcoeffs u
      vcoeffs = tr (L.inv gvandermonde) L.#> evalcoeffs v
   in Grid $
        fromList
          [ (ucoeffs `atIndex` i) * (vcoeffs `atIndex` j)
            | i <- [0 .. np - 1],
              j <- [0 .. np - 1]
          ]

--------------------------------------------------------------------------------

-- | Derivative of a basis function
deriv1 :: (Eq a, Num a, Storable a, U.Unbox a) => Int -> Vector a
deriv1 i = poly2vector $ deriv $ monomial (fromIntegral i) 1

-- | Derivatives of all basis functions
derivs1 :: (Element a, Eq a, Num a, U.Unbox a) => Matrix a
derivs1 = fromColumns [deriv1 i | i <- [0 .. np -1]]

derivs :: (Eq a, Field a, U.Unbox a) => Matrix a
derivs = gvandermonde <> derivs1 <> L.inv gvandermonde

-- derivs2 :: (Eq a, Field a, U.Unbox a) => Matrix a
-- derivs2 = derivs <> derivs

expfilter1 :: (Container Vector a, RealFloat a) => Matrix a
expfilter1 =
  diag $
    fromList
      [ let n0 = np `div` 2
         in if i <= n0
              then 1
              else peps ** (fromIntegral (i - n0) / fromIntegral (np - n0))
        | i <- [0 .. np -1]
      ]

expfilters :: (Container Vector a, Field a, RealFloat a, U.Unbox a) => Matrix a
expfilters = gvandermonde <> expfilter1 <> L.inv gvandermonde

bnd :: (Element a, Fractional a) => Matrix a
bnd =
  diag $
    fromList
      [ if i == 0
          then -1 / fromIntegral (np - 1)
          else if i == np - 1 then 1 / fromIntegral (np - 1) else 0
        | i <- [0 .. np - 1]
      ]

weights :: (Element a, Fractional a) => Matrix a
weights =
  diag $
    fromList
      [ if i == 0 || i == np - 1
          then 1 / fromIntegral (np - 1)
          else 2 / fromIntegral (np - 1)
        | i <- [0 .. np - 1]
      ]

--------------------------------------------------------------------------------

mkOp :: (Element a, Num a) => ((Int, Int) -> a) -> ((Int, Int) -> a) -> Op a
mkOp fx fy =
  Op $
    fromLists
      [ [ fx (ir, ic) * fy (jr, jc)
          | ic <- [0 .. np - 1],
            jc <- [0 .. np - 1]
        ]
        | ir <- [0 .. np - 1],
          jr <- [0 .. np - 1]
      ]

vandermonde :: (Container Vector a, Eq a, Floating a, U.Unbox a) => Op a
vandermonde = mkOp (gvandermonde `atIndex`) (gvandermonde `atIndex`)

expfilter :: (Container Vector a, Field a, RealFloat a, U.Unbox a) => Op a
expfilter = mkOp (expfilters `atIndex`) (expfilters `atIndex`)

delta :: Field a => Op a
delta = mkOp delta1 delta1
  where
    delta1 (i, j) = if i == j then 1 else 0

diagonal :: (Element a, Num a) => Grid a -> Op a
diagonal (Grid xs) = Op $ diag xs

originR :: Field a => Op a
originR =
  Op $
    fromLists
      [ [ if ic == ir && jc == jr && ic == jc then 1 else 0
          | ic <- [0 .. np - 1],
            jc <- [0 .. np - 1]
        ]
        | ir <- [0 .. np - 1],
          jr <- [0 .. np - 1]
      ]

valueAtOriginR :: (Eq a, Field a, U.Unbox a) => Op a
valueAtOriginR =
  Op $
    fromRows
      [ let t = coordT ! (ir, jr)
            r = 0
            v = (t + r) / 2
            u = (t - r) / 2
         in getGrid (gevalCoeffs (u, v))
        | ir <- [0 .. np - 1],
          jr <- [0 .. np - 1]
      ]

reverseU :: (Element a, Num a) => Op a
reverseU = mkOp reverse1 delta1
  where
    reverse1 (i, j) = if i == np-1-j then 1 else 0
    delta1 (i, j) = if i == j then 1 else 0

reverseV :: (Element a, Num a) => Op a
reverseV = mkOp delta1 reverse1
  where
    reverse1 (i, j) = if i == np-1-j then 1 else 0
    delta1 (i, j) = if i == j then 1 else 0

-- change (t,r) -> (t,-r); (u,v) -> (v,u)
reverseR :: (Element a, Num a) => Op a
reverseR =
  Op $
    fromLists
      [ [ if ic == jr && jc == ir then 1 else 0
          | ic <- [0 .. np - 1],
            jc <- [0 .. np - 1]
        ]
        | ir <- [0 .. np - 1],
          jr <- [0 .. np - 1]
      ]

symmetrizeU :: Field a => Op a
symmetrizeU = (1 / 2) *^ (delta ^+^ reverseU)

symmetrizeV :: Field a => Op a
symmetrizeV = (1 / 2) *^ (delta ^+^ reverseV)

symmetrizeR :: Field a => Op a
symmetrizeR = (1 / 2) *^ (delta ^+^ reverseR)

-- | Derivatives of a grid
derivU :: (Eq a, Field a, U.Unbox a) => Op a
derivU = mkOp (derivs `atIndex`) delta1
  where
    delta1 (i, j) = if i == j then 1 else 0

derivV :: (Eq a, Field a, U.Unbox a) => Op a
derivV = mkOp delta1 (derivs `atIndex`)
  where
    delta1 (i, j) = if i == j then 1 else 0

derivUU :: (Eq a, Field a, U.Unbox a) => Op a
derivUU = derivU * derivU

derivUV :: (Eq a, Field a, U.Unbox a) => Op a
derivUV = derivU * derivV

derivVV :: (Eq a, Field a, U.Unbox a) => Op a
derivVV = derivV * derivV

-- | Combined derivatives
-- > t = v + u
-- > r = v - u
-- > t^2 = v^2 + 2 v u + u^2
-- > r^2 = v^2 - 2 v u + u^2
derivT :: (Eq a, Field a, U.Unbox a) => Op a
derivT = (derivV ^+^ derivU) ^/ 2

derivR :: (Eq a, Field a, U.Unbox a) => Op a
derivR = (derivV ^-^ derivU) ^/ 2

derivTT :: (Eq a, Field a, U.Unbox a) => Op a
derivTT = derivT * derivT

derivRR :: (Eq a, Field a, U.Unbox a) => Op a
derivRR = derivR * derivR

-- | 1/r dr
r1derivR :: (Field a, RealFloat a, U.Unbox a) => Op a
r1derivR =
  ((delta - originR) * diagonal coordR1 * derivR + originR * derivRR)
    * symmetrizeR

-- t = v + u   v = (t + r) / 2
-- r = v - u   u = (t - r) / 2

-- | 1/r du
r1derivU :: (Field a, RealFloat a, U.Unbox a) => Op a
r1derivU =
  ((delta - originR) * diagonal coordR1 * derivU - originR * derivUU)
    * symmetrizeR

-- | 1/r dv
r1derivV :: (Field a, RealFloat a, U.Unbox a) => Op a
r1derivV =
  ((delta - originR) * diagonal coordR1 * derivV + originR * derivVV)
    * symmetrizeR

-- | 1/r^2
r2project :: (Field a, RealFloat a, U.Unbox a) => Op a
r2project =
  ( (delta - originR) * diagonal (gmap (^ 2) coordR1)
      + originR * ((1 / 2) *^ derivRR)
  )
    * (delta - valueAtOriginR)
    * symmetrizeR

-- | zero at origin
zproject :: (Eq a, Field a, U.Unbox a) => Op a
zproject = (delta - valueAtOriginR) * symmetrizeR

-- | Boundary locations
bndU :: (Container Vector a, Fractional a) => Op a
bndU = mkOp (bnd `atIndex`) delta1
  where
    delta1 (i, j) = if i == j then 1 else 0

bndV :: (Container Vector a, Fractional a) => Op a
bndV = mkOp delta1 (bnd `atIndex`)
  where
    delta1 (i, j) = if i == j then 1 else 0

metric :: (Container Vector a, Fractional a) => Op a
metric = mkOp (weights `atIndex`) (weights `atIndex`)

--------------------------------------------------------------------------------

solve :: Field a => Op a -> Grid a -> Grid a
solve (Op a) (Grid x) =
  Grid $ head $ toColumns $ fromJust $ linearSolve a (asColumn x)

solveLS :: Field a => Op a -> Grid a -> Grid a
solveLS (Op a) (Grid x) =
  Grid $ head $ toColumns $ linearSolveLS a (asColumn x)

--------------------------------------------------------------------------------

newtype MaxAbs a = MaxAbsInternal a
  deriving (Eq, Ord, Read, Show)

maxabs :: Num a => a -> MaxAbs a
maxabs x = MaxAbsInternal (abs x)

getMaxAbs :: MaxAbs a -> a
getMaxAbs (MaxAbsInternal x) = x

instance Ord a => Semigroup (MaxAbs a) where
  MaxAbsInternal x <> MaxAbsInternal y = MaxAbsInternal $ max x y

instance (Num a, Ord a) => Monoid (MaxAbs a) where
  mempty = MaxAbsInternal 0
