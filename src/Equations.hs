{-# LANGUAGE BangPatterns #-}
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

eqnG01 :: Num a => State a -> State a -> State a -> State a -> State a -> State a -> a
{-# INLINE eqnG01 #-}
eqnG01 !s@(State g q w) !r2s@(State r2g r2q _) !r1sr@(State r1gr r1qr r1wr) !su@(State gu qu wu) !sv@(State gv qv wv) !suv@(State guv quv wuv) =
  -2 * r2q + 2 * r2g - 2 * r1qr + quv

eqnG22 :: Num a => State a -> State a -> State a -> State a -> State a -> State a -> a
{-# INLINE eqnG22 #-}
eqnG22 !s@(State g q w) !r2s@(State r2g r2q _) !r1sr@(State r1gr r1qr r1wr) !su@(State gu qu wu) !sv@(State gv qv wv) !suv@(State guv quv wuv) =
  2 * q * g ^ 2 * r1qr + g ^ 2 * qv * qu - 4 * q ^ 2 * g ^ 2 * wv * wu + 2 * q ^ 2 * gv * gu - 2 * q * g ^ 2 * quv - 2 * q ^ 2 * g * guv

eqnSW2 :: Num a => State a -> State a -> State a -> State a -> State a -> State a -> a
{-# INLINE eqnSW2 #-}
eqnSW2 !s@(State g q w) !r2s@(State r2g r2q _) !r1sr@(State r1gr r1qr r1wr) !su@(State gu qu wu) !sv@(State gv qv wv) !suv@(State guv quv wuv) =
  4 * q * r1wr - 2 * wv * qu - 2 * qv * wu - 4 * q * wuv

-- Constraints

eqnG00 :: Num a => Coord a -> State a -> State a -> State a -> State a -> State a -> a
{-# INLINE eqnG00 #-}
eqnG00 !c@(Coord u v) !s@(State g q w) !r2s@(State r2g r2q r2w) !r1su@(State r1gu r1qu r1wu) !su@(State gu qu wu) !suu@(State guu quu wuu) =
  let r = v - u
   in 4 * q ^ 2 * r1qu - 4 * q ^ 2 * r1gu + 4 * q * r1qu * g - 4 * q * r1gu * g - 4 * r * q * r2q * qu + 4 * r * q * r2g * qu + g * qu ^ 2 - 4 * q ^ 2 * g * wu ^ 2 - 4 * r * q * r2q * gu + 4 * r * q * r2g * gu + 2 * q * qu * gu - 2 * q * g * quu

eqnG11 :: Num a => Coord a -> State a -> State a -> State a -> State a -> State a -> a
{-# INLINE eqnG11 #-}
eqnG11 !c@(Coord u v) !s@(State g q w) !r2s@(State r2g r2q r2w) !r1sv@(State r1gv r1qv r1wv) !sv@(State gv qv wv) !svv@(State gvv qvv wvv) =
  let r = v - u
   in -4 * q ^ 2 * r1qv + 4 * q ^ 2 * r1gv - 4 * q * r1qv * g + 4 * q * r1gv * g + 4 * r * q * r2q * qv - 4 * r * q * r2g * qv + g * qv ^ 2 - 4 * q ^ 2 * g * wv ^ 2 + 4 * r * q * r2q * gv - 4 * r * q * r2g * gv + 2 * q * qv * gv - 2 * q * g * qvv
