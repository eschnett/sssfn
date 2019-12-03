{-# OPTIONS_GHC -fno-warn-type-defaults #-}

-- Ideas taken from <http://sriku.org/blog/2019/03/12/automatic-differentiation-dual-numbers--taylor-numbers/>

module Dual
  ( Dual (..),
  )
where

import Control.Applicative
import Data.VectorSpace
import Foreign
import GHC.Generics
import qualified Test.QuickCheck as QC

default (Int)

data Dual a = Dual {primal :: a, dual :: a}
  deriving (Read, Show, Generic, Foldable, Functor, Traversable)

instance (Bounded a, Num a) => Bounded (Dual a) where

  minBound = Dual minBound 0

  maxBound = Dual maxBound 0

instance QC.Arbitrary a => QC.Arbitrary (Dual a) where

  arbitrary = QC.applyArbitrary2 Dual

  shrink = QC.genericShrink

instance Storable a => Storable (Dual a) where

  sizeOf _ = 2 * sizeOf (undefined :: a)

  alignment _ = alignment (undefined :: a)

  peek ptr = do
    x <- peekElemOff (castPtr ptr) 0
    dx <- peekElemOff (castPtr ptr) 1
    return (Dual x dx)

  poke ptr (Dual x dx) = do
    pokeElemOff (castPtr ptr) 0 x
    pokeElemOff (castPtr ptr) 1 dx

instance Applicative Dual where

  pure x = Dual x x

  liftA2 f (Dual x dx) (Dual y dy) = Dual (f x y) (f dx dy)

instance Num a => Num (Dual a) where

  (+) = liftA2 (+)

  negate = fmap negate

  Dual x dx * Dual y dy = Dual (x * y) (x * dy + dx * y)

  -- signum = fmap signum
  signum (Dual x dx) = Dual (signum x) 0

  -- abs = fmap abs
  abs (Dual x dx) = Dual (abs x) (dx * signum x)

  fromInteger i = Dual (fromInteger i) 0

instance Eq a => Eq (Dual a) where
  -- Equality ignores the dual part
  Dual x _ == Dual y _ = x == y

instance Ord a => Ord (Dual a) where
  -- Ordering ignores the dual part
  compare (Dual x _) (Dual y _) = compare x y

instance Fractional a => Fractional (Dual a) where

  fromRational r = Dual (fromRational r) 0

  Dual x dx / Dual y dy = Dual (x / y) (dx / y - (x / y) * (dy / y))

instance Floating a => Floating (Dual a) where

  pi = Dual pi 0

  exp (Dual x dx) = Dual (exp x) (dx * exp x)

  log (Dual x dx) = Dual (log x) (dx / x)

  sqrt (Dual x dx) = Dual (sqrt x) (dx / (2 * sqrt x))

  Dual x dx ** Dual y dy = Dual (x ** y) (dx * (y / x) + dy * log x)

  -- logBase = _

  sin (Dual x dx) = Dual (sin x) (dx * cos x)

  cos (Dual x dx) = Dual (cos x) (- dx * sin x)

  tan (Dual x dx) = Dual (tan x) (dx / cos x ^ 2)

  asin (Dual x dx) = Dual (asin x) (dx / sqrt (1 - x ^ 2))

  acos (Dual x dx) = Dual (acos x) (- dx / sqrt (1 - x ^ 2))

  atan (Dual x dx) = Dual (atan x) (dx / sqrt (1 + x ^ 2))

  sinh (Dual x dx) = Dual (sinh x) (dx * cosh x)

  cosh (Dual x dx) = Dual (cosh x) (dx * sinh x)

  tanh (Dual x dx) = Dual (tanh x) (dx / cosh x ^ 2)

  asinh (Dual x dx) = Dual (asinh x) (dx / sqrt (1 + x ^ 2))

  acosh (Dual x dx) = Dual (acosh x) (dx / sqrt (x ^ 2 -1))

  atanh (Dual x dx) = Dual (atanh x) (dx / sqrt (1 - x ^ 2))

-- log1mexp = _
-- log1p = _
-- expm1 = _
-- log1pexp = _
-- log1mexp = _

instance AdditiveGroup a => AdditiveGroup (Dual a) where

  (^+^) = liftA2 (^+^)

  negateV = fmap negateV

  zeroV = pure zeroV

instance VectorSpace a => VectorSpace (Dual a) where

  type Scalar (Dual a) = Scalar a

  (*^) a = fmap (a *^)
