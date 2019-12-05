{-# OPTIONS_GHC -fno-warn-type-defaults #-}

module Equations where

import Control.Applicative
import Data.VectorSpace
import Foreign
import GHC.Generics
import qualified Test.QuickCheck as QC

default (Int)

--------------------------------------------------------------------------------

-- | Variables
data Var = VarG | VarQ | VarW
  deriving (Bounded, Enum, Eq, Ord, Read, Show)

-- | Equations
data Eqn = EqnG01 | EqnG22 | EqnSW2 | EqnG00 | EqnG11
  deriving (Bounded, Enum, Eq, Ord, Read, Show)

--------------------------------------------------------------------------------

data Coord a = Coord {u :: !a, v :: !a}
  deriving (Eq, Ord, Read, Show, Generic, Foldable, Functor, Traversable)

instance Storable a => Storable (Coord a) where

  sizeOf _ = 2 * sizeOf (undefined :: a)

  alignment _ = alignment (undefined :: a)

  peek ptr = do
    u <- peekElemOff (castPtr ptr) 0
    v <- peekElemOff (castPtr ptr) 1
    return (Coord u v)

  poke ptr (Coord u v) = do
    pokeElemOff (castPtr ptr) 0 u
    pokeElemOff (castPtr ptr) 1 v

instance QC.Arbitrary a => QC.Arbitrary (Coord a) where

  arbitrary = QC.applyArbitrary2 Coord

  shrink = QC.genericShrink

instance Applicative Coord where

  pure x = Coord x x

  liftA2 f (Coord u1 v1) (Coord u2 v2) = Coord (f u1 u2) (f v1 v2)

--------------------------------------------------------------------------------

data State a = State {g :: !a, q :: !a, w :: !a}
  deriving (Eq, Ord, Read, Show, Generic, Foldable, Functor, Traversable)

instance Storable a => Storable (State a) where

  sizeOf _ = 3 * sizeOf (undefined :: a)

  alignment _ = alignment (undefined :: a)

  peek ptr = do
    g <- peekElemOff (castPtr ptr) 0
    q <- peekElemOff (castPtr ptr) 1
    w <- peekElemOff (castPtr ptr) 2
    return (State g q w)

  poke ptr (State g q w) = do
    pokeElemOff (castPtr ptr) 0 g
    pokeElemOff (castPtr ptr) 1 q
    pokeElemOff (castPtr ptr) 2 w

instance QC.Arbitrary a => QC.Arbitrary (State a) where

  arbitrary = QC.applyArbitrary3 State

  shrink = QC.genericShrink

instance Applicative State where

  pure x = State x x x

  liftA2 f (State g1 q1 w1) (State g2 q2 w2) =
    State (f g1 g2) (f q1 q2) (f w1 w2)

instance Num a => Num (State a) where

  (+) = liftA2 (+)

  negate = fmap negate

  (*) = liftA2 (*)

  abs = fmap abs

  signum = fmap signum

  fromInteger i = pure (fromInteger i)

instance Num a => AdditiveGroup (State a) where

  (^+^) = (+)

  zeroV = pure 0

  negateV = negate

instance Num a => VectorSpace (State a) where

  type Scalar (State a) = a

  (*^) a = fmap (a *)

--------------------------------------------------------------------------------

-- Evolution equations

eqnG01 :: Fractional a => Coord a -> State a -> State a -> State a -> State a -> a
{-# INLINE eqnG01 #-}
eqnG01 (Coord u v) (State g q w) (State gu qu wu) (State gv qv wv) (State guv quv wuv) = (-14 * q * (u - v) + 4 * (u - v) ^ 2 * qv - 2 * gv - 4 * u ^ 2 * qu + 8 * u * v * qu - 4 * v ^ 2 * qu + 2 * gu + u ^ 3 * quv - 3 * u ^ 2 * v * quv + 3 * u * v ^ 2 * quv - v ^ 3 * quv - u * guv + v * guv) / ((u - v) * (q * (u - v) ^ 2 - g))

eqnG22 :: Fractional a => Coord a -> State a -> State a -> State a -> State a -> a
{-# INLINE eqnG22 #-}
eqnG22 (Coord u v) (State g q w) (State gu qu wu) (State gv qv wv) (State guv quv wuv) = (4 * q ^ 4 * (u - v) ^ 7 * (-1 + (u - v) ^ 2 * wv * wu) + 2 * q * (u - v) * g * (4 * g ^ 2 + (u - v) ^ 2 * gv * (3 * (u - v) ^ 2 * qu + gu) + (u - v) ^ 4 * qv * ((u - v) ^ 2 * qu + 3 * gu) - 2 * (u - v) ^ 3 * g * (4 * qv - 4 * qu + (u - v) * quv)) + 4 * q ^ 3 * (u - v) ^ 5 * (-4 * g + (u - v) * ((u - v) ^ 2 * qv - gv - u ^ 2 * qu + 2 * u * v * qu - v ^ 2 * qu + gu + u ^ 3 * quv - 3 * u ^ 2 * v * quv + 3 * u * v ^ 2 * quv - v ^ 3 * quv)) - q ^ 2 * (u - v) ^ 3 * (4 * g ^ 2 * (-7 + 2 * (u - v) ^ 2 * wv * wu) + (u - v) ^ 4 * qv * (3 * (u - v) ^ 2 * qu + gu) + (u - v) ^ 2 * gv * ((u - v) ^ 2 * qu + 3 * gu) - 2 * (u - v) * g * (3 * (u - v) ^ 2 * qv + 5 * gv - 3 * u ^ 2 * qu + 6 * u * v * qu - 3 * v ^ 2 * qu - 5 * gu - 2 * u * guv + 2 * v * guv)) + g ^ 2 * (4 * (u - v) * g ^ 2 * wv * wu - (u - v) * ((u - v) ^ 2 * qv * (3 * (u - v) ^ 2 * qu + gu) + gv * ((u - v) ^ 2 * qu + 3 * gu)) - 2 * g * ((u - v) ^ 2 * qv - gv - u ^ 2 * qu + 2 * u * v * qu - v ^ 2 * qu + gu - 2 * u * guv + 2 * v * guv))) / (4 * (u - v) * (q * (u - v) ^ 2 - g) * (q * (u - v) ^ 2 + g) ^ 3)

eqnSW2 :: Fractional a => Coord a -> State a -> State a -> State a -> State a -> a
{-# INLINE eqnSW2 #-}
eqnSW2 (Coord u v) (State g q w) (State gu qu wu) (State gv qv wv) (State guv quv wuv) = - ((u - v) * (((u - v) ^ 2 * qv - gv) * wu + wv * ((u - v) ^ 2 * qu - gu)) + 2 * q * (u - v) ^ 2 * (2 * wv - 2 * wu + (u - v) * wuv) - 2 * g * (wv - wu + (u - v) * wuv)) / (2 * (u - v) * (q * (u - v) ^ 2 - g) * (q * (u - v) ^ 2 + g))

-- Constraints

eqnG00 :: Fractional a => Coord a -> State a -> State a -> State a -> a
{-# INLINE eqnG00 #-}
eqnG00 (Coord u v) (State g q w) (State gu qu wu) (State guu quu wuu) = - (4 * q ^ 3 * (u - v) ^ 4 * (-2 + (u - v) ^ 2 * wu ^ 2) - q * (4 * g ^ 2 * (5 + (u - v) ^ 2 * wu ^ 2) - 12 * (u - v) * g * ((u - v) ^ 2 * qu + gu) + (u - v) ^ 2 * (3 * (u - v) ^ 4 * qu ^ 2 - 2 * (u - v) ^ 2 * qu * gu - gu ^ 2)) + 2 * q ^ 2 * (u - v) ^ 2 * (-2 * g * (-5 + (u - v) ^ 2 * wu ^ 2) + (u - v) * (-2 * (u - v) ^ 2 * qu - 2 * gu + (u - v) * ((u - v) ^ 2 * quu - guu))) + g * ((u - v) ^ 4 * qu ^ 2 + 4 * g ^ 2 * wu ^ 2 + 2 * (u - v) ^ 2 * qu * gu - 3 * gu ^ 2 - 2 * g * (8 * (u - v) * qu + (u - v) ^ 2 * quu - guu))) / (2 * (- (q * (u - v) ^ 2) + g) ^ 2 * (q * (u - v) ^ 2 + g))

eqnG11 :: Fractional a => Coord a -> State a -> State a -> State a -> a
{-# INLINE eqnG11 #-}
eqnG11 (Coord u v) (State g q w) (State gv qv wv) (State gvv qvv wvv) = - (4 * q ^ 3 * (u - v) ^ 4 * (-2 + (u - v) ^ 2 * wv ^ 2) - q * (4 * g ^ 2 * (5 + (u - v) ^ 2 * wv ^ 2) + 12 * (u - v) * g * ((u - v) ^ 2 * qv + gv) + (u - v) ^ 2 * (3 * (u - v) ^ 4 * qv ^ 2 - 2 * (u - v) ^ 2 * qv * gv - gv ^ 2)) + 2 * q ^ 2 * (u - v) ^ 2 * (-2 * g * (-5 + (u - v) ^ 2 * wv ^ 2) + (u - v) * (2 * (u - v) ^ 2 * qv + 2 * gv + (u - v) * ((u - v) ^ 2 * qvv - gvv))) + g * ((u - v) ^ 4 * qv ^ 2 + 4 * g ^ 2 * wv ^ 2 + 2 * (u - v) ^ 2 * qv * gv - 3 * gv ^ 2 - 2 * g * (-8 * (u - v) * qv + (u - v) ^ 2 * qvv - gvv))) / (2 * (- (q * (u - v) ^ 2) + g) ^ 2 * (q * (u - v) ^ 2 + g))

-- Evolution equations on axis

eqnAxisG01 :: Fractional a => Coord a -> State a -> State a -> State a -> State a -> a
{-# INLINE eqnAxisG01 #-}
eqnAxisG01 (Coord u v) (State g q w) (State gu qu wu) (State gv qv wv) (State guv quv wuv) = (-14 * q + (u - v) ^ 2 * quv - guv) / (q * (u - v) ^ 2 - g)

eqnAxisG22 :: Fractional a => Coord a -> State a -> State a -> State a -> State a -> a
{-# INLINE eqnAxisG22 #-}
eqnAxisG22 (Coord u v) (State g q w) (State gu qu wu) (State gv qv wv) (State guv quv wuv) = (4 * q ^ 4 * (u - v) ^ 6 * (-4 + (u - v) ^ 2 * wv ^ 2 + 2 * (u - v) ^ 2 * wv * wu + (u - v) ^ 2 * wu ^ 2) + 16 * q ^ 3 * (u - v) ^ 4 * (-4 * g + (u - v) ^ 4 * quv) + 2 * q * g * (16 * g ^ 2 + (u - v) ^ 2 * ((u - v) ^ 4 * qv ^ 2 + gv ^ 2 + u ^ 4 * qu ^ 2 - 4 * u ^ 3 * v * qu ^ 2 + 6 * u ^ 2 * v ^ 2 * qu ^ 2 - 4 * u * v ^ 3 * qu ^ 2 + v ^ 4 * qu ^ 2 + 6 * u ^ 2 * qu * gu - 12 * u * v * qu * gu + 6 * v ^ 2 * qu * gu + gu ^ 2 + 2 * gv * (3 * (u - v) ^ 2 * qu + gu) + 2 * (u - v) ^ 2 * qv * (3 * gv + (u - v) ^ 2 * qu + 3 * gu)) - 8 * (u - v) ^ 4 * g * quv) + g ^ 2 * (-3 * (u - v) ^ 4 * qv ^ 2 - 3 * gv ^ 2 - 2 * u ^ 2 * gv * qu + 4 * u * v * gv * qu - 2 * v ^ 2 * gv * qu - 3 * u ^ 4 * qu ^ 2 + 12 * u ^ 3 * v * qu ^ 2 - 18 * u ^ 2 * v ^ 2 * qu ^ 2 + 12 * u * v ^ 3 * qu ^ 2 - 3 * v ^ 4 * qu ^ 2 + 4 * g ^ 2 * (wv + wu) ^ 2 - 6 * gv * gu - 2 * u ^ 2 * qu * gu + 4 * u * v * qu * gu - 2 * v ^ 2 * qu * gu - 3 * gu ^ 2 - 2 * (u - v) ^ 2 * qv * (gv + 3 * (u - v) ^ 2 * qu + gu) + 16 * g * guv) - q ^ 2 * (8 * (u - v) ^ 2 * g ^ 2 * (-14 + (u - v) ^ 2 * wv ^ 2 + 2 * (u - v) ^ 2 * wv * wu + (u - v) ^ 2 * wu ^ 2) + (u - v) ^ 4 * (3 * (u - v) ^ 4 * qv ^ 2 + 3 * gv ^ 2 + 3 * u ^ 4 * qu ^ 2 - 12 * u ^ 3 * v * qu ^ 2 + 18 * u ^ 2 * v ^ 2 * qu ^ 2 - 12 * u * v ^ 3 * qu ^ 2 + 3 * v ^ 4 * qu ^ 2 + 2 * u ^ 2 * qu * gu - 4 * u * v * qu * gu + 2 * v ^ 2 * qu * gu + 3 * gu ^ 2 + 2 * (u - v) ^ 2 * qv * (gv + 3 * (u - v) ^ 2 * qu + gu) + 2 * gv * ((u - v) ^ 2 * qu + 3 * gu)) + 16 * (u - v) ^ 4 * g * guv)) / (16 * (q * (u - v) ^ 2 - g) * (q * (u - v) ^ 2 + g) ^ 3)

eqnAxisSW2 :: Fractional a => Coord a -> State a -> State a -> State a -> State a -> a
{-# INLINE eqnAxisSW2 #-}
eqnAxisSW2 (Coord u v) (State g q w) (State gu qu wu) (State gv qv wv) (State guv quv wuv) = (gv * wu - u ^ 2 * qu * wu + 2 * u * v * qu * wu - v ^ 2 * qu * wu - (u - v) ^ 2 * qv * (wv + wu) + wu * gu + wv * (gv - (u - v) ^ 2 * qu + gu) - 4 * q * u ^ 2 * wuv + 8 * q * u * v * wuv - 4 * q * v ^ 2 * wuv + 4 * g * wuv) / (4 * (q * (u - v) ^ 2 - g) * (q * (u - v) ^ 2 + g))

-- Constraints on axis

eqnAxisG00 :: Fractional a => Coord a -> State a -> State a -> State a -> State a -> a
{-# INLINE eqnAxisG00 #-}
eqnAxisG00 (Coord u v) (State g q w) (State gu qu wu) (State gv qv wv) (State guu quu wuu) = - (4 * q ^ 3 * (u - v) ^ 4 * (-8 + (u - v) ^ 2 * wv ^ 2 + 2 * (u - v) ^ 2 * wv * wu + (u - v) ^ 2 * wu ^ 2) - q * (4 * g ^ 2 * (20 + (u - v) ^ 2 * wv ^ 2 + 2 * (u - v) ^ 2 * wv * wu + (u - v) ^ 2 * wu ^ 2) - 24 * (u - v) * g * ((u - v) ^ 2 * qv + gv + u ^ 2 * qu - 2 * u * v * qu + v ^ 2 * qu + gu) + (u - v) ^ 2 * (3 * (u - v) ^ 4 * qv ^ 2 - gv ^ 2 + 3 * u ^ 4 * qu ^ 2 - 12 * u ^ 3 * v * qu ^ 2 + 18 * u ^ 2 * v ^ 2 * qu ^ 2 - 12 * u * v ^ 3 * qu ^ 2 + 3 * v ^ 4 * qu ^ 2 + 2 * (u - v) ^ 2 * qv * (- gv + 3 * (u - v) ^ 2 * qu - gu) - 2 * u ^ 2 * qu * gu + 4 * u * v * qu * gu - 2 * v ^ 2 * qu * gu - gu ^ 2 - 2 * gv * ((u - v) ^ 2 * qu + gu))) + g * ((u - v) ^ 4 * qv ^ 2 - 3 * gv ^ 2 + 2 * u ^ 2 * gv * qu - 4 * u * v * gv * qu + 2 * v ^ 2 * gv * qu + u ^ 4 * qu ^ 2 - 4 * u ^ 3 * v * qu ^ 2 + 6 * u ^ 2 * v ^ 2 * qu ^ 2 - 4 * u * v ^ 3 * qu ^ 2 + v ^ 4 * qu ^ 2 + 4 * g ^ 2 * (wv + wu) ^ 2 - 6 * gv * gu + 2 * u ^ 2 * qu * gu - 4 * u * v * qu * gu + 2 * v ^ 2 * qu * gu - 3 * gu ^ 2 + 2 * (u - v) ^ 2 * qv * (gv + (u - v) ^ 2 * qu + gu) - 8 * g * (4 * (u - v) * qv + 4 * (u - v) * qu + u ^ 2 * quu - 2 * u * v * quu + v ^ 2 * quu - guu)) + 4 * q ^ 2 * (u - v) ^ 2 * (- (g * (-20 + (u - v) ^ 2 * wv ^ 2 + 2 * (u - v) ^ 2 * wv * wu + (u - v) ^ 2 * wu ^ 2)) + 2 * (u - v) * (- ((u - v) ^ 2 * qv) - gv - u ^ 2 * qu + 2 * u * v * qu - v ^ 2 * qu - gu + u ^ 3 * quu - 3 * u ^ 2 * v * quu + 3 * u * v ^ 2 * quu - v ^ 3 * quu - u * guu + v * guu))) / (8 * (- (q * (u - v) ^ 2) + g) ^ 2 * (q * (u - v) ^ 2 + g))

eqnAxisG11 :: Fractional a => Coord a -> State a -> State a -> State a -> State a -> a
{-# INLINE eqnAxisG11 #-}
eqnAxisG11 (Coord u v) (State g q w) (State gu qu wu) (State gv qv wv) (State gvv qvv wvv) = - (4 * q ^ 3 * (u - v) ^ 4 * (-8 + (u - v) ^ 2 * wv ^ 2 + 2 * (u - v) ^ 2 * wv * wu + (u - v) ^ 2 * wu ^ 2) + g * ((u - v) ^ 4 * qv ^ 2 - 3 * gv ^ 2 + 2 * u ^ 2 * gv * qu - 4 * u * v * gv * qu + 2 * v ^ 2 * gv * qu + u ^ 4 * qu ^ 2 - 4 * u ^ 3 * v * qu ^ 2 + 6 * u ^ 2 * v ^ 2 * qu ^ 2 - 4 * u * v ^ 3 * qu ^ 2 + v ^ 4 * qu ^ 2 - 8 * g * (-4 * (u - v) * qv + (u - v) ^ 2 * qvv - gvv - 4 * u * qu + 4 * v * qu) + 4 * g ^ 2 * (wv + wu) ^ 2 - 6 * gv * gu + 2 * u ^ 2 * qu * gu - 4 * u * v * qu * gu + 2 * v ^ 2 * qu * gu - 3 * gu ^ 2 + 2 * (u - v) ^ 2 * qv * (gv + (u - v) ^ 2 * qu + gu)) + 4 * q ^ 2 * (u - v) ^ 2 * (- (g * (-20 + (u - v) ^ 2 * wv ^ 2 + 2 * (u - v) ^ 2 * wv * wu + (u - v) ^ 2 * wu ^ 2)) + 2 * (u - v) * ((u - v) ^ 2 * qv + gv + u ^ 3 * qvv - 3 * u ^ 2 * v * qvv + 3 * u * v ^ 2 * qvv - v ^ 3 * qvv - u * gvv + v * gvv + u ^ 2 * qu - 2 * u * v * qu + v ^ 2 * qu + gu)) - q * (4 * g ^ 2 * (20 + (u - v) ^ 2 * wv ^ 2 + 2 * (u - v) ^ 2 * wv * wu + (u - v) ^ 2 * wu ^ 2) + 24 * (u - v) * g * ((u - v) ^ 2 * qv + gv + u ^ 2 * qu - 2 * u * v * qu + v ^ 2 * qu + gu) + (u - v) ^ 2 * (3 * (u - v) ^ 4 * qv ^ 2 - gv ^ 2 + 3 * u ^ 4 * qu ^ 2 - 12 * u ^ 3 * v * qu ^ 2 + 18 * u ^ 2 * v ^ 2 * qu ^ 2 - 12 * u * v ^ 3 * qu ^ 2 + 3 * v ^ 4 * qu ^ 2 + 2 * (u - v) ^ 2 * qv * (- gv + 3 * (u - v) ^ 2 * qu - gu) - 2 * u ^ 2 * qu * gu + 4 * u * v * qu * gu - 2 * v ^ 2 * qu * gu - gu ^ 2 - 2 * gv * ((u - v) ^ 2 * qu + gu)))) / (8 * (- (q * (u - v) ^ 2) + g) ^ 2 * (q * (u - v) ^ 2 + g))
