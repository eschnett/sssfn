{-# LANGUAGE TypeApplications #-}

import SSSFN

main :: IO ()
main =
  do
    putStrLn "Spherically Symmetric Scalar Field in double-Null coordinates"
    putStrLn $ show (gcoordU @Double)
    putStrLn $ show (derivU @Double)
    putStrLn $ show (derivU #> gcoordU @Double)
    putStrLn "Done."
  where
    ini (u, v) = cos (pi * u)
