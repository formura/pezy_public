{-# LANGUAGE FlexibleContexts #-}

import Text.Regex.TDFA
import Data.List
import Control.Monad
import Text.Printf

calc start end lines = do
    lines <- return $ dropWhile (not . (=~ start)) lines
    lines <- return $ takeWhile (not . (=~ end)) lines

    -- print lines

    let ops = map words $ filter (=~ "^\t[a-zA-Z0-9]+\\.[a-zA-Z0-9]+") lines

    let total = length ops

    ops <- return $ sortOn (negate . length) $ groupBy (\a b -> head a == head b) $ sortOn head ops

    forM_ ops $ \op -> do
        when (length op >= 10) $ printf "%-12s: %d\n" (head $ head op) (length op)
    putStrLn ""

    -- print ops
    -- print $ length ops

    let madOps    = ["d.mad", "d.nmad", "d.msub"]
    let singleOps = ["d.add", "d.sub", "d.mul"]

    let madOpNum    = sum [ length ls | op <- madOps, let Just ls = find (\s -> head (head s) == op) ops ]
    let singleOpNum = sum [ length ls | op <- singleOps, let Just ls = find (\s -> head (head s) == op) ops ]

    printf   "mad ops:   %d\n" madOpNum
    printf   "sngle ops: %d\n" singleOpNum
    printf   "dp total:  %d = %d + %d\n" (madOpNum + singleOpNum) madOpNum singleOpNum
    printf   "flop:      %d = %d * 2 + %d\n" (madOpNum * 2 + singleOpNum) madOpNum singleOpNum
    printf   "other:     %d\n" (total - (madOpNum + singleOpNum))
    printf   "total:     %d\n" total

    return ()

main = do
    lines <- lines <$> getContents

    printf "Precalc:\n"
    calc "// precalc start" "// precalc end" lines

    printf "\n"

    printf "Step:\n"
    calc "// step start" "// step end" lines
