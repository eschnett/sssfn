{-# OPTIONS_GHC -fno-warn-type-defaults #-}

-- Ideas taken from <http://sriku.org/blog/2019/03/12/automatic-differentiation-dual-numbers--taylor-numbers/>

module DualF
  ( DualF (..),
    Dual,
    mkDual,
    dual,
  )
where

import Control.Applicative
import Data.Functor.Identity
import Data.VectorSpace
import Foreign
import GHC.Generics
import qualified Test.QuickCheck as QC

default (Int)

data DualF v a = Dual {primal :: a, dualf :: v a}
  deriving (Read, Show, Generic, Foldable, Functor, Traversable)

instance (Bounded a, AdditiveGroup (v a)) => Bounded (DualF v a) where

  minBound = Dual minBound zeroV

  maxBound = Dual maxBound zeroV

instance (QC.Arbitrary a, QC.Arbitrary (v a)) => QC.Arbitrary (DualF v a) where

  arbitrary = QC.applyArbitrary2 Dual

  shrink = QC.genericShrink

instance (Storable a, Storable (v a)) => Storable (DualF v a) where

  sizeOf _ = sizeOf (undefined :: a) + sizeOf (undefined :: v a)

  alignment _ = alignment (undefined :: a) `max` alignment (undefined :: v a)

  peek ptr = do
    x <- peek (castPtr ptr)
    dx <- peekByteOff (castPtr ptr) (sizeOf x)
    return (Dual x dx)

  poke ptr (Dual x dx) = do
    poke (castPtr ptr) x
    pokeByteOff (castPtr ptr) (sizeOf x) dx

instance Applicative v => Applicative (DualF v) where

  pure x = Dual x (pure x)

  liftA2 f (Dual x dx) (Dual y dy) = Dual (f x y) (liftA2 f dx dy)

instance (Num a, VectorSpace (v a), Scalar (v a) ~ a) => Num (DualF v a) where

  Dual x dx + Dual y dy = Dual (x + y) (dx ^+^ dy)

  negate (Dual x dx) = Dual (negate x) (negateV dx)

  Dual x dx * Dual y dy = Dual (x * y) (x *^ dy ^+^ dx ^* y)

  signum (Dual x dx) = Dual (signum x) zeroV

  abs (Dual x dx) = Dual (abs x) (signum x *^ dx)

  fromInteger i = Dual (fromInteger i) zeroV

instance Eq a => Eq (DualF v a) where
  -- Equality ignores the dual part
  Dual x _ == Dual y _ = x == y

instance Ord a => Ord (DualF v a) where
  -- Ordering ignores the dual part
  compare (Dual x _) (Dual y _) = compare x y

instance
  (Fractional a, VectorSpace (v a), Scalar (v a) ~ a) =>
  Fractional (DualF v a)
  where

  fromRational r = Dual (fromRational r) zeroV

  Dual x dx / Dual y dy = Dual (x / y) (dx ^/ y ^-^ (x / y ^ 2) *^ dy)

instance
  (Floating a, VectorSpace (v a), Scalar (v a) ~ a) =>
  Floating (DualF v a)
  where

  pi = Dual pi zeroV

  exp (Dual x dx) = Dual (exp x) (dx ^* exp x)

  log (Dual x dx) = Dual (log x) (dx ^/ x)

  sqrt (Dual x dx) = Dual (sqrt x) (dx ^/ (2 * sqrt x))

  Dual x dx ** Dual y dy = Dual (x ** y) (dx ^* (y / x) ^+^ dy ^* log x)

  -- logBase = _

  sin (Dual x dx) = Dual (sin x) (dx ^* cos x)

  cos (Dual x dx) = Dual (cos x) (dx ^* (- sin x))

  tan (Dual x dx) = Dual (tan x) (dx ^/ cos x ^ 2)

  asin (Dual x dx) = Dual (asin x) (dx ^/ sqrt (1 - x ^ 2))

  acos (Dual x dx) = Dual (acos x) (dx ^/ (- sqrt (1 - x ^ 2)))

  atan (Dual x dx) = Dual (atan x) (dx ^/ sqrt (1 + x ^ 2))

  sinh (Dual x dx) = Dual (sinh x) (dx ^* cosh x)

  cosh (Dual x dx) = Dual (cosh x) (dx ^* sinh x)

  tanh (Dual x dx) = Dual (tanh x) (dx ^/ cosh x ^ 2)

  asinh (Dual x dx) = Dual (asinh x) (dx ^/ sqrt (1 + x ^ 2))

  acosh (Dual x dx) = Dual (acosh x) (dx ^/ sqrt (x ^ 2 -1))

  atanh (Dual x dx) = Dual (atanh x) (dx ^/ sqrt (1 - x ^ 2))

-- log1mexp = _
-- log1p = _
-- expm1 = _
-- log1pexp = _
-- log1mexp = _

instance
  (AdditiveGroup a, AdditiveGroup (v a)) =>
  AdditiveGroup (DualF v a)
  where

  Dual x dx ^+^ Dual y dy = Dual (x ^+^ y) (dx ^+^ dy)

  negateV (Dual x dx) = Dual (negateV x) (negateV dx)

  zeroV = Dual zeroV zeroV

instance
  (VectorSpace a, VectorSpace (v a), Scalar (v a) ~ Scalar a) =>
  VectorSpace (DualF v a)
  where

  type Scalar (DualF v a) = Scalar a

  a *^ Dual x dx = Dual (a *^ x) (a *^ dx)

--------------------------------------------------------------------------------

type Dual = DualF Identity

mkDual :: a -> a -> Dual a
mkDual x dx = Dual x (Identity dx)

dual :: Dual a -> a
dual = runIdentity . dualf
