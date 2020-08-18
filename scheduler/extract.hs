#!/usr/bin/env stack script --resolver=lts-11.0

{-# LANGUAGE ViewPatterns #-}

import Text.Regex.TDFA
import Control.Monad
import Text.Printf
import Data.Char
import Data.List.Split (splitOn)
import Data.Maybe
import Data.List
import Text.ParserCombinators.Parsec hiding (State)
import qualified Text.Parsec.Token as P
import Text.Parsec.Language (javaStyle)
import Data.Either
import qualified Data.Map as Map
import Data.Map.Merge.Lazy
import Control.Monad.State
import qualified Data.Set as Set
import Control.Exception.Base (assert)
import Data.Unique

-- (+1 * q[i - 2][j + 0][k + 0] - 8 * q[i - 1][j + 0][k + 0] + 8 * q[i + 1][j + 0][k + 0] - 1 * q[i + 2][j + 0][k + 0]) * r12h;

data DExpr = DExpr
    { stencil :: [(Int, (Int, Int, Int))]
    , factor :: (Int, Int)
    }
    deriving (Show)

removeSpaces :: String -> String
removeSpaces (a: b: c: cs)
    | isSpace b && not (isAlpha a && isAlpha c) =
        removeSpaces (a: c: cs)
removeSpaces (c:cs) =
    c: removeSpaces cs
removeSpaces [] = []

parseFunc :: String -> (String, DExpr)
parseFunc s =
    let re_fn = "double (d_.+)\\(.*\\){return(\\(.+\\))/\\(([0-9]+)((\\*h)+)\\);}"
        re_body = "([+-]?[0-9]+)\\*q\\[i([+-][0-9]+)\\]\\[j([+-][0-9]+)\\]\\[k([+-][0-9]+)\\]"
        [ma] = s =~ re_fn :: [[String]]
        body = ma !! 2
        expr = body =~ re_body :: [[String]]
        f "" = 1
        f cs = read . g . filter (not . isSpace) $ cs
        g ('+': cs) = cs
        g cs = cs
    in (ma!!1,
        DExpr
            { stencil = map (\ss -> (f $ ss!!1, (f $ ss!!2, f $ ss!!3, f $ ss!!4))) expr
            , factor = (read $ ma!!3, length (ma!!4) `div` 2)
            })

parseUse :: [(String, DExpr)] -> String -> (String, [(String, (Int, Int, Int))])
parseUse funcs s =
    let re_use = "double (.+) = (d_.+)\\(ix, iy, iz, (.+)+\\)"
        [ma] = s =~ re_use
        var = ma!!1
        f = ma!!2
        buf = ma!!3
        refs = stencil $ fromJust (lookup f funcs)
    in (var, [(buf, ref) | (_, ref) <- refs])

lexer = P.makeTokenParser javaStyle

lexeme = P.lexeme lexer
whiteSpace = P.whiteSpace lexer
symbol = P.symbol lexer
naturalOrFloat = P.naturalOrFloat lexer
identifier = P.identifier lexer
parens = P.parens lexer
braces = P.braces lexer
brackets = P.brackets lexer
commaSep = P.commaSep lexer
commaSep1 = P.commaSep1 lexer

type Ident = String

data Function
    = Function Ident [Ident] Expr
    deriving (Show)

data Statement
    = Decl (Maybe Ident) Expr
    deriving (Show, Eq, Ord)

data Ref
    = ArrayRef Ident [Int] [Expr]
    | VarRef Ident
    deriving (Show, Eq, Ord)

data Expr
    = BinOp Op Expr Expr
    | UniOp Op Expr
    | Call Ident [Expr]
    | Const Value
    | VarExpr Ref
    | AssignExpr Expr Expr
    deriving (Show, Eq, Ord)

data Op = Add | Sub | Mul | Div
    deriving (Show, Eq, Ord)

data Value = VInt Int | VDouble Double
    deriving (Show, Eq, Ord)

ppStmt :: Statement -> String
ppStmt stmt = case stmt of
    Decl (Just var) rhs -> "let " ++ var ++ " = " ++ ppExpr rhs ++ ";"
    Decl Nothing rhs -> ppExpr rhs ++ ";"

ppExpr :: Expr -> String
ppExpr expr = case expr of
    BinOp op l r -> "(" ++ ppExpr l ++ " " ++ ppOp op ++ " " ++ ppExpr r ++ ")"
    UniOp op e -> ppOp op ++ " " ++ ppExpr e
    Call f args -> f ++ "(" ++ intercalate ", " (map ppExpr args) ++ ")"
    Const (VInt n) -> show n
    Const (VDouble d) -> show d
    VarExpr ref -> ppRef ref

ppRef :: Ref -> String
ppRef ref = case ref of
    ArrayRef v _dims ixs -> v ++ concat [ "[" ++ ppExpr ix ++ "]" | ix <- ixs ]
    VarRef v -> v

ppOp :: Op -> String
ppOp Add = "+"
ppOp Sub = "-"
ppOp Mul = "*"
ppOp Div = "/"

opname :: Op -> String
opname Add = "add"
opname Sub = "sub"
opname Mul = "mul"
opname Div = "div"

mapStmt :: (Expr -> Expr) -> Statement -> Statement
mapStmt f stmt = case stmt of
    Decl v e -> Decl v $ f e

getExpr :: Statement -> Expr
getExpr stmt = case stmt of
    Decl _ e -> e

-----

parseFunction :: Parser Function
parseFunction = do
    whiteSpace
    void $ symbol "double"
    name <- identifier
    args <- parens $ commaSep $ do
        _ty <- identifier
        var <- identifier
        _ixs <- many $ brackets parseExpr
        return var

    body <- braces $ do
        void $ symbol "return"
        e <- parseExpr
        void $ symbol ";"
        return e

    return $ Function name args body

-- TODO: parse from text
dims = reverse [c_NX+2*c_Ns, c_NY+2*c_Ns, c_NZ+2*c_Ns]

parseAssign :: Parser Statement
parseAssign = whiteSpace >> (declStmt <|> assignStmt) where
    declStmt = do
        void $ symbol "double"
        ident <- identifier
        symbol "="
        rhs <- parseExpr
        symbol ";"
        return $ Decl (Just ident) rhs

    assignStmt = do
        ident <- identifier
        ixs <- many $ brackets parseExpr
        symbol "="
        rhs <- parseExpr
        symbol ";"
        return $ Decl Nothing $ AssignExpr (VarExpr $ ArrayRef ident dims ixs) rhs

parseExpr :: Parser Expr
parseExpr = expr where
    expr = term `chainl1` addop
    term = uniop `chainl1` mulop

    uniop = do
        f <- option id $ (UniOp Sub <$ symbol "-") <|> (id <$ symbol "+")
        e <- fact
        return $ f e

    fact = const <|> try call <|> var <|> parens expr

    const = do
        v <- naturalOrFloat
        return $ Const $ either (VInt . fromIntegral) VDouble v

    call = do
        ident <- identifier
        args <- parens $ commaSep1 $ expr
        return $ Call ident args

    var = do
        v <- identifier
        ixs <- many $ brackets expr
        return $ VarExpr $ if null ixs then VarRef v else ArrayRef v dims ixs

    addop = BinOp <$> (Add <$ symbol "+" <|> Sub <$ symbol "-")
    mulop = BinOp <$> (Mul <$ symbol "*" <|> Div <$ symbol "/")

constProp :: Expr -> Expr
constProp e = case e of
    BinOp op (constProp -> l) (constProp -> r) -> case (l, r) of
        (Const cl, Const cr) -> Const $ evalBin op cl cr
        _ -> BinOp op l r

    UniOp op (constProp -> e) -> case e of
        Const ce -> Const $ evalUni op ce
        _ -> UniOp op e

    Call f (map constProp -> ixs) ->
        -- TODO: Should be eval for const?
        Call f ixs

    Const _ -> e

    VarExpr (VarRef var) -> case lookup var constantTable of
        Just val -> Const val
        _ -> e

    VarExpr (ArrayRef var dims ixs) ->
        VarExpr $ ArrayRef var dims $ map constProp ixs

    AssignExpr refe e ->
        AssignExpr (constProp refe) (constProp e)

toDouble :: Value -> Double
toDouble (VInt n) = fromIntegral n
toDouble (VDouble d) = d

evalBin :: Op -> Value -> Value -> Value
evalBin op l r = case (op, l, r) of
    (Add, VInt a, VInt b) -> VInt $ a + b
    (Add, a, b) -> VDouble $ toDouble a + toDouble b

    (Sub, VInt a, VInt b) -> VInt $ a - b
    (Sub, a, b) -> VDouble $ toDouble a - toDouble b

    (Mul, VInt a, VInt b) -> VInt $ a + b
    (Mul, a, b) -> VDouble $ toDouble a * toDouble b

    (Div, a, b) -> VDouble $ toDouble a / toDouble b

evalUni :: Op -> Value -> Value
evalUni op e = case (op, e) of
    (Sub, VInt a) -> VInt $ -a
    (Sub, VDouble d) -> VDouble $ -d

weaken :: Expr -> Expr
weaken e = case e of
    BinOp op (weaken -> l) (weaken -> r) -> case (op, l, r) of
        (Div, _, Const v) -> BinOp Mul l $ Const $ VDouble $ 1.0 / toDouble v

        (Add, _, Const (VInt 0)) -> l
        (Add, _, Const (VDouble d)) | abs d < 1e-10 -> l
        (Add, Const (VInt 0), _) -> r
        (Add, Const (VDouble d), _) | abs d < 1e-10 -> r

        (Mul, _, Const (VInt 1)) -> l
        (Mul, _, Const (VDouble d)) | abs (d - 1) < 1e-10 -> l
        (Mul, Const (VInt 1), _) -> r
        (Mul, Const (VDouble d), _) | abs (d - 1) < 1e-10 -> r

        _ -> BinOp op l r

    UniOp op (weaken -> e) ->
        UniOp op e

    Call f (map weaken -> ixs) -> case (f, ixs) of
        ("pow", [a, Const (VDouble n)]) ->
            if n > 0
                then foldl1 (BinOp Mul) $ replicate (round n) a
                else foldl1 (BinOp Mul) $ replicate (round $ abs n) $ BinOp Div (Const $ VDouble 1.0) a
        ("pow", ixs) -> error $ show ixs
        _ -> Call f ixs

    Const _ -> e

    VarExpr v -> case v of
        ArrayRef v dims ixs ->
            let wixs = map weaken ixs
                e = snd $
                        foldl1 (\(dim, wix) (_, e) -> (0, BinOp Add (BinOp Mul wix $ Const $ VInt dim) e)) $
                        zip (tail $ scanr (*) 1 dims) wixs
                vr = VarExpr (VarRef v)
            in BinOp Add vr e

        _ -> e

    AssignExpr refe e ->
        AssignExpr (weaken refe) (weaken e)

subExprs :: Expr -> Map.Map Expr Int
subExprs e = case e of
    BinOp op l r ->
        Map.singleton e 1 `mergeSubs` subExprs l `mergeSubs` subExprs r

    UniOp op e ->
        Map.singleton e 1 `mergeSubs` subExprs e

    Call f ixs ->
        Map.singleton e 1 `mergeSubs` (foldr mergeSubs Map.empty $ map subExprs ixs)

    Const _ -> Map.singleton e 1
    VarExpr _ -> Map.singleton e 1

removeDiv :: Expr -> Expr
removeDiv e = case e of
    BinOp op (removeDiv -> l) (removeDiv -> r) -> case (op, l, r) of
        (Div, Const (VDouble 1.0), VarExpr (VarRef "r")) ->
            VarExpr (VarRef "r_recip")
        (Div, _, _) ->
            error "Unexpected div operator"

        _ ->
            BinOp op l r

    UniOp op (removeDiv -> l) ->
        UniOp op l

    Call f ixs ->
        Call f (map removeDiv ixs)

    Const _ -> e
    VarExpr _ -> e

mergeSubs :: Map.Map Expr Int -> Map.Map Expr Int -> Map.Map Expr Int
mergeSubs = merge preserveMissing preserveMissing (zipWithMatched $ \_ x y -> x + y)

inlining :: [(Ident, Function)] -> Expr -> Expr
inlining dict e = go e where
    go e = case e of
        Call f ixs ->
            let Just (Function _ args body) = lookup f dict in
            subst (zip args ixs) body

        BinOp op (go -> l) (go -> r) ->
            BinOp op l r
        UniOp op (go -> l) ->
            UniOp op l
        Const _ -> e
        VarExpr _ -> e
        AssignExpr (go -> var) (go -> r) ->
            AssignExpr var r

subst :: [(Ident, Expr)] -> Expr -> Expr
subst dict e = go e where
    go e = case e of
        BinOp op (go -> l) (go -> r) ->
            BinOp op l r
        UniOp op (go -> l) ->
            UniOp op l
        Call f ixs ->
            Call f $ map go ixs
        Const _ -> e

        VarExpr (VarRef v) ->
            case lookup v dict of
                Just w -> w
                Nothing -> e

        VarExpr (ArrayRef v dims ixs) ->
            let Just (VarExpr (VarRef w)) = lookup v dict in
            VarExpr $ ArrayRef w dims $ map go ixs

type LProgram = [(Int, LInst)]

-- TODO: Support immediate argument
data LInst
    = LOp LOp [Int]
    | LInt Int
    | LDouble Double
    | LVar Ident
    | LLoad Int
    | LStore Int Int
    deriving (Show)

data LOp
    = IAdd
    | ISub
    | IMul
    | IDiv
    | DAdd
    | DSub
    | DMul
    | DDiv
    | DNeg
    deriving (Show)

isIntOp :: LOp -> Bool
isIntOp op = case op of
    IAdd -> True
    ISub -> True
    IMul -> True
    IDiv -> True
    DAdd -> True
    DSub -> True
    DMul -> True
    DDiv -> True
    DNeg -> True

-- FIXME:
opToLOp :: Op -> LOp
opToLOp op = case op of
    Add -> IAdd
    Sub -> ISub
    Mul -> IMul
    Div -> IDiv

ppInst :: Int -> LInst -> String
ppInst ix inst = case inst of
    LOp op args ->
        printf "{\"type\": \"alu\", \"loc\": %d, \"opc\": %s, \"args\": %s}" ix (show op) (show args)
    LInt n ->
        printf "{\"type\": \"i32\", \"loc\": %d, \"value\": %d}" ix n
    -- TODO: encode to hex
    LDouble d ->
        printf "{\"type\": \"f64\", \"loc\": %d, \"value\": %f}" ix d
    LVar v ->
        printf "{\"type\": \"var\", \"loc\": %d, \"ident\": %s}" ix (show v)
    LLoad addr ->
        printf "{\"type\": \"load\", \"loc\": %d, \"addr\": %d}" ix addr
    LStore addr e ->
        printf "{\"type\": \"store\", \"addr\": %d, \"addr\": %d}" addr e

instRefs :: LInst -> [Int]
instRefs inst = case inst of
    LOp _ args -> args
    LInt _ -> []
    LDouble _ -> []
    LVar _ -> []
    LLoad addr -> [addr]
    LStore addr val -> [addr, val]

type Lower = State LowerState

data LowerState
    = LowerState
        { lsVars :: Map.Map Ident Int
        , lsInsts :: [(Int, LInst)] -- Instructions (reverse order)
        , lsCnt :: Int -- Cache inst size
        }

lowerExpr :: Expr -> Lower Int
lowerExpr e = case e of
    BinOp op l r -> do
        ll <- lowerExpr l
        lr <- lowerExpr r
        emitInst $ LOp (opToLOp op) [ll, lr]

    UniOp op e -> do
        le <- lowerExpr e
        emitInst $ LOp (opToLOp op) [le]

    Const (VInt n) ->
        emitInst $ LInt $ fromIntegral n

    Const (VDouble d) ->
        emitInst $ LDouble d

    VarExpr (VarRef v) -> do
        vars <- gets lsVars
        case Map.lookup v vars of
            Nothing ->
                emitInst $ LVar v
            Just ix ->
                return ix

    VarExpr _ ->
        undefined

    AssignExpr var e -> do
        lvar <- lowerExpr var
        le <- lowerExpr e
        emitInst $ LStore lvar le

lowerStmt :: Statement -> Lower ()
lowerStmt stmt = case stmt of
    Decl (Just v) e -> do
        le <- lowerExpr e
        ls <- get
        let vars = lsVars ls
        let nvars = Map.insert v le vars
        put $ ls { lsVars = nvars }

    Decl Nothing e -> do
        void $ lowerExpr e

findVar :: Ident -> Lower (Maybe Int)
findVar e = do
    vars <- gets lsVars
    return $ Map.lookup e vars

emitInst :: LInst -> Lower Int
emitInst inst = do
    ls <- get
    let insts = lsInsts ls
    let ix = lsCnt ls
    put $ ls { lsInsts = (ix, inst) : insts, lsCnt = ix + 1 }
    return ix

doLower :: [Statement] -> LProgram
doLower stmts =
    let initState = LowerState Map.empty [] 0
    in reverse $ lsInsts $ execState (mapM_ lowerStmt stmts) initState

type CodeGen = StateT CGState IO

data CGState
    = CGState
        { cgLocation :: Map.Map Int Loc
        , cgFreeRegs :: [Reg]
        , cgUsedRegs :: [Reg]
        , cgRestRefs :: Map.Map Int Int
        , cgMC :: [String]
        }

data Loc
    = LocStack Int
    | LocReg Reg

data Reg
    = Reg RegType String
    deriving (Eq)

ppReg :: Reg -> String
ppReg (Reg _ name) = name

data RegType
    = RTInt
    | RTDouble
    deriving (Eq)

registers :: [Reg]
registers =
    [ Reg RTInt (show name) | name <- [1..30] ] ++
    [ Reg RTDouble (show name) | name <- [0..31] ]

zr :: Reg
zr = Reg RTInt "r0"

emitMC :: String -> CodeGen ()
emitMC mc = do
    modify $ \state -> state { cgMC = mc : cgMC state }

cgLoad :: Int -> CodeGen Reg
cgLoad ix = do
    st <- get
    ret <- case Map.lookup ix (cgLocation st) of
        Nothing -> undefined

        Just (LocReg r) -> do
            cgTouch r
            return r

        Just (LocStack stk) -> do
            r <- cgAlloc RTInt -- TODO: Fix
            emitMC $ printf "// TODO: load from stack"
            return r

    return ret

cgTouch :: Reg -> CodeGen ()
cgTouch reg = do
    st <- get
    let regs = cgUsedRegs st
    let nregs = reg : (regs \\ [reg])
    put $ st { cgUsedRegs = nregs }

cgAlloc :: RegType -> CodeGen Reg
cgAlloc = undefined

cgLoadVar :: Ident -> CodeGen Reg
cgLoadVar = undefined

cgInst :: LInst -> CodeGen (Maybe Reg)
cgInst inst = case inst of
    LOp op args -> do
        largs <- mapM cgLoad args
        if isIntOp op
            then do
            dest <- cgAlloc RTInt
            case (op, largs) of
                (IAdd, [a1, a2]) -> do
                    emitMC $ printf "i.add %s %s %s" (ppReg dest) (ppReg a1) (ppReg a2)
                (ISub, [a1, a2]) -> do
                    emitMC $ printf "i.sub %s %s %s" (ppReg dest) (ppReg a1) (ppReg a2)
                (IMul, [a1, a2]) -> do
                    emitMC $ printf "i.mul %s %s %s" (ppReg dest) (ppReg a1) (ppReg a2)
                (IDiv, [a1, a2]) -> do
                    emitMC $ printf "i.div %s %s %s" (ppReg dest) (ppReg a1) (ppReg a2)
                _ -> undefined
            return $ Just dest

            else do
            dest <- cgAlloc RTDouble
            case (op, largs) of
                (DAdd, [a1, a2]) -> do
                    emitMC $ printf "d.add %s %s %s" (ppReg dest) (ppReg a1) (ppReg a2)
                (DSub, [a1, a2]) -> do
                    emitMC $ printf "d.sub %s %s %s" (ppReg dest) (ppReg a1) (ppReg a2)
                (DMul, [a1, a2]) -> do
                    emitMC $ printf "d.mul %s %s %s" (ppReg dest) (ppReg a1) (ppReg a2)
                (DDiv, [a1, a2]) -> do
                    emitMC $ printf "d.div %s %s %s" (ppReg dest) (ppReg a1) (ppReg a2)
                _ -> undefined
            return $ Just dest

    LInt v -> do
        dest <- cgAlloc RTDouble
        if v >= 0x10000000
            then do
            emitMC $ printf "i.movhi %s %d" (ppReg dest) (v `div` 0x10000)
            emitMC $ printf "i.ori %s %s %d" (ppReg dest) (ppReg dest) (v `mod` 0x10000)
            else do
            emitMC $ printf "i.ori %s %s %d" (ppReg dest) (ppReg zr) (0x10000 :: Int)
        return $ Just dest

    LDouble v -> do
        dest <- cgAlloc RTDouble
        emitMC $ printf "d.movi %f; TODO" v
        return $ Just dest

    LVar var -> do
        var <- cgLoadVar var
        dest <- cgAlloc RTInt
        emitMC $ printf "i.ldd %s %s" (ppReg dest) (ppReg var)
        return $ Just dest

    LLoad addr -> do
        arg <- cgLoad addr
        dest <- cgAlloc RTDouble
        emitMC $ printf "d.eldd %s %s 0" (ppReg dest) (ppReg arg)
        return $ Just dest

    LStore addr val -> do
        lval <- cgLoad val
        laddr <- cgLoad addr
        emitMC $ printf "d.esd %s 0 %s" (ppReg laddr) (ppReg lval)
        return Nothing

codeGen :: [(Int, LInst)] -> IO [String]
codeGen insts = do
    let initState = CGState
            { cgLocation = Map.empty
            , cgFreeRegs = registers
            , cgUsedRegs = []
            , cgRestRefs = Map.empty
            , cgMC = []
            }
    st <- execStateT (mapM_ cgInst $ map snd insts) initState
    return $ reverse $ cgMC st

-----

c_NX = 10
c_NY = 10
c_NZ = 10
c_NT = 5
c_Ns = 2
c_MX = 10
c_MY = 10
c_MZ = 10
c_T_MAX = 100

c_h = 1.0 / (fromIntegral c_NX * fromIntegral c_MX)
c_cfl = 1.0
c_dt = c_cfl * c_h * c_h
c_re = 150.0
c_mu = 1.0
c_gm = 1.4
c_c = c_mu / c_re
c_c2 = (c_gm - 1) * c_mu / c_re

constantTable =
    [ ("NX", VInt c_NX)
    , ("NY", VInt c_NY)
    , ("NZ", VInt c_NZ)
    , ("NT", VInt c_NT)
    , ("Ns", VInt c_Ns)

    , ("MX", VInt c_MX)
    , ("MY", VInt c_MY)
    , ("MZ", VInt c_MZ)

    , ("T_MAX", VInt c_T_MAX)

    , ("h", VDouble c_h)
    , ("cfl", VDouble c_cfl)
    , ("dt", VDouble c_dt)
    , ("re", VDouble c_re)
    , ("mu", VDouble c_mu)
    , ("gm", VDouble c_gm)
    , ("c", VDouble c_c)
    , ("c2", VDouble c_c2)
    ]

main :: IO ()
main = do
    ls <- lines <$> getContents
    ls <- return $ filter (not . null) ls
    let [fs, us, es] = splitOn ["-----"] ls

    -- putStrLn $ show $ map removeSpaces fs

    let funcs = map (parseFunc . removeSpaces) fs
    -- forM_ funcs print

    -- let funcs = map (either (\e -> error $ "parse failed: " ++ show e) id . parse parseFunction "") fs
    -- forM_ funcs print

    let pu :: String -> (String, String, String)
        pu s =
            let re_use = "double (.+)=(d_.+)\\(ix,iy,iz,(.)p_buf\\);"
                [ma] = s =~ re_use
            in (ma!!1, ma!!2, ma!!3)

    let uses = map (pu . removeSpaces) us
    -- forM_ uses print

    let grp = groupBy (\a b -> fst a == fst b) $
                sortOn (\(s, _) -> (length s, s)) $
                map (\(v, f, e) -> (f, e)) uses

    forM_ grp $ \g -> do
        putStrLn $ "// " ++ fst (head g) ++ ": " ++ show (map snd g)
    putStrLn ""

    let opnums =
            [ ((func, ss, fields), (length ss - 1 + sum [ if abs fa == 1 then 0 else 1 | (fa, _) <- ss ] + 1) * length fields)
            | g <- grp
            , let func = fst (head g)
            , let fields = map snd g
            , let DExpr { stencil = ss, factor = (n, h) } = fromJust $ lookup func funcs
            ]

    printf "// Total ops: %d\n" $ sum $ map snd opnums
    -- forM_ opnums $ \opn ->
    --     printf "// - %s\n" $ show opn
    putStrLn ""

    -- forM_ [(x, y, z) | x <- [-2..2::Int], y <- [-2 .. -2::Int], z <- [-2 .. -2::Int]] $ \(x, y, z) -> do
    --     printf "const auto *e%d%d%d = &buf.p[ix%+d][iy%+d][iz%+d];\n" (x+2) (y+2) (z+2) x y z
    -- putStrLn ""

    -- putStrLn "if (get_pid() == 2048) { ix--; }"
    -- putStrLn ""

    -- forM_ grp $ \g -> do
    --     let func = fst (head g)
    --     let fields = map snd g

    --     let code = printf "const auto e_%s = (%s) * consts.rev%dh%d;"
    --                     (drop 2 func) body n h
    --         DExpr { stencil = ss, factor = (n, h) } = fromJust $ lookup func funcs

    --         -- body = concat
    --         --         [ case abs f of
    --         --             1  -> printf "%c %s" sig ff :: String
    --         --             -- 2  -> printf "%c (e%d%d%d + e%d%d%d)" sig (x+2) (y+2) (z+2) (x+2) (y+2) (z+2)
    --         --             _ -> printf "%c consts.num%d * %s" sig (abs f) ff
    --         --         | (f, (x, y, z)) <- ss
    --         --         , let sig = if f > 0 then '+' else '-'
    --         --         , let ff = printf "e%d00[%d*(NZ+2*Ns)+%d]" (x+2) (y+2) (z+2) :: String
    --         --         ]

    --         body = concat
    --                 [ case abs f of
    --                     1  -> printf "+ (%s)" body :: String
    --                     _ -> printf " + consts.num%d * (%s)" (abs f) body
    --                 | g@((f, _):_) <-
    --                     groupBy (\a b -> abs (fst a) == abs (fst b)) $
    --                     sortOn (abs . fst) $  ss
    --                 , let body = concatMap (\(f, (x, y, z)) ->
    --                         let ff = printf "e%d00[%d*(NZ+2*Ns)+%d]" (x+2) (y+2) (z+2) :: String
    --                         in printf " %c %s" (if f > 0 then '+' else '-') ff :: String) g
    --                 ]

    --     putStrLn $ code

    --     forM_ fields $ \[field] -> do
    --         -- printf "const double %c_%s = e_%s.%c;\n" field (drop 2 func) (drop 2 func) field
    --         printf "#define %c_%s e_%s.%c\n" field (drop 2 func) (drop 2 func) field


    ----------
    -- フィールド内回し

    forM_ [(x, y, z) | x <- [-2..2::Int], y <- [-2 .. -2::Int], z <- [-2 .. -2::Int]] $ \(x, y, z) -> do
        printf "const double *e%d%d%d = (double*)&buf.p[ix%+d][iy%+d][iz%+d];\n" (x+2) (y+2) (z+2) x y z
    putStrLn ""

    forM_ grp $ \g -> do
        let func = fst (head g)
        let fields = map snd g

        forM_ fields $ \[field] -> do
            let code = printf "const double %c_%s = (%s) * consts.rev%dh%d;"
                            field (drop 2 func) (dropWhile (== '+') body) n h
                DExpr { stencil = ss, factor = (n, h) } = fromJust $ lookup func funcs

                Just fieldIx = field `lookup` (zip "ruvwp" [0 :: Int ..])

                body = concat
                        [ case f of
                            1  -> printf "+%s" ff :: String
                            -1 -> printf "-%s" ff
                            _ -> printf "%cconsts.num%d*%s" (if f > 0 then '+' else '-') (abs f) ff
                        | (f, (x, y, z)) <- ss
                        , let ff = printf "e%d00[(%d*(NZ+2*Ns)+%d)*5+%d]" (x+2) (y+2) (z+2) fieldIx :: String
                        ]

                -- body = concat
                --         [ case abs f of
                --             1  -> printf "+ (%s)" body :: String
                --             _ -> printf "+ consts.num%d * (%s)" (abs f) body
                --         | g@((f, _):_) <-
                --             groupBy (\a b -> abs (fst a) == abs (fst b)) $
                --             sortOn (abs . fst) $  ss

                --         , let body' = concatMap (\(f, (x, y, z)) ->
                --                 let ff = printf "e%d00[(%d*(NZ+2*Ns)+%d)*5 + %d]" (x+2) (y+2) (z+2) fieldIx :: String
                --                 in printf "%c %s" (if f > 0 then '+' else '-') ff :: String) g
                --         , let body = if head body' == '+' then tail body' else body'
                --         ]

            putStrLn code

        putStrLn ""
        putStrLn "if (get_pid() == 2048) { return; }"
        putStrLn ""

    ----------
    -- まとめて生成

    -- forM_ [(x, y, z) | x <- [-2..2::Int], y <- [-2 .. -2::Int], z <- [-2 .. -2::Int]] $ \(x, y, z) -> do
    --     printf "const auto *e%d%d%d = (double)&buf.p[ix%+d][iy%+d][iz%+d];\n" (x+2) (y+2) (z+2) x y z
    -- putStrLn ""

    -- forM_ "pruvw" $ \field -> do
    --     forM_ grp $ \g -> do
    --         let func = fst (head g)
    --         let fields = map snd g
    --         when ([field] `elem` fields) $ do

    --             let code = printf "const double %c_%s = (%s) * consts.rev%dh%d;"
    --                             field (drop 2 func) (dropWhile (== '+') body) n h
    --                 DExpr { stencil = ss, factor = (n, h) } = fromJust $ lookup func funcs

    --                 Just fieldIx = field `lookup` (zip "ruvwp" [0 :: Int ..])

    --                 body = concat
    --                         [ case f of
    --                             1  -> printf "+%s" ff :: String
    --                             -1 -> printf "-%s" ff
    --                             _ -> printf "%cconsts.num%d*%s" (if f > 0 then '+' else '-') (abs f) ff
    --                         | (f, (x, y, z)) <- ss
    --                         , let ff = printf "e%d00[(%d*(NZ+2*Ns)+%d)*5+%d]" (x+2) (y+2) (z+2) fieldIx :: String
    --                         ]

    --                 -- body = concat
    --                 --         [ case abs f of
    --                 --             1  -> printf "+ (%s)" body :: String
    --                 --             _ -> printf "+ consts.num%d * (%s)" (abs f) body
    --                 --         | g@((f, _):_) <-
    --                 --             groupBy (\a b -> abs (fst a) == abs (fst b)) $
    --                 --             sortOn (abs . fst) $  ss

    --                 --         , let body' = concatMap (\(f, (x, y, z)) ->
    --                 --                 let ff = printf "e%d00[(%d*(NZ+2*Ns)+%d)].%c" (x+2) (y+2) (z+2) field :: String
    --                 --                 in printf "%c %s" (if f > 0 then '+' else '-') ff :: String) g
    --                 --         , let body = if head body' == '+' then tail body' else body'
    --                 --         ]
    --             putStrLn code

    --     putStrLn ""
    --     putStrLn "if (get_pid() == 2048) { ix = 0; }"
    --     putStrLn ""

    ----------
    -- アセンブリ生成

    -- forM_ [(x, y, z) | x <- [-2..2::Int], y <- [-2 .. -2::Int], z <- [-2 .. -2::Int]] $ \(x, y, z) -> do
    --     printf "const double *e%d%d%d = (double*)&buf.p[ix%+d][iy%+d][iz%+d];\n" (x+2) (y+2) (z+2) x y z
    -- putStrLn ""

    -- forM_ grp $ \g -> do
    --     let func = fst (head g)
    --     let fields = map snd g
    --     let Just (DExpr { stencil = ss, factor = (n, h) }) = lookup func funcs
    --     let grpStencil =
    --             map (\g -> (abs $ fst $ head g, g)) $
    --             groupBy (\a b -> abs (fst a) == abs (fst b)) $
    --             sortOn (\(f, _) -> abs f) ss

    --     let genVar = do
    --             u <- newUnique
    --             let ix = hashUnique u
    --             let var = printf "t%d" ix :: String
    --             printf "double %s;\n" var
    --             return var

    --     let nz = 36
    --     let ns = 2

    --     forM_ fields $ \[field] -> do
    --         let resultVar = printf "%c_%s" field (drop 2 func) :: String
    --         printf "// Code for %s\n" resultVar

    --         let Just fi = field `lookup` (zip "ruvwp" [0 :: Int ..])

    --         let cell x y z =
    --                 let Just fi = field `lookup` (zip "ruvwp" [0 :: Int ..])
    --                 in printf "e%d00[(%d*(NZ+2*Ns)+%d) * 5 + %d]" (x+2) (y+2) (z+2) fi :: String

    --         terms <- forM grpStencil $ \(sfact, (f, (x, y, z)): poss) -> do
    --             var <- genVar
    --             if f > 0
    --                 then printf "%s = %s;\n" var (cell x y z)
    --                 else printf "%s = -%s;\n" var (cell x y z)
    --             let go prev poss = case poss of
    --                     [] -> return prev
    --                     (f, (x, y, z)): rest -> do
    --                         -- nv <- genVar
    --                         -- printf "d_eldd(%s, e%d00, ((%d*(NZ+2*Ns)+%d) * 5 + %d) * 8);\n" nv (x+2) (y+2) (z+2) fi
    --                         nv <- (return (cell x y z :: String) :: IO String)
    --                         if f > 0
    --                             then printf "d_add(%s, %s, %s);\n" prev prev nv
    --                             else printf "d_sub(%s, %s, %s);\n" prev prev nv
    --                         go prev rest
    --             res <- go var poss
    --             printf "d_mul(%s, %s, consts.num%d);\n" res res sfact :: IO ()
    --             return res

    --         let go (Just res) [] = return res
    --             go Nothing (v:vs) = go (Just v) vs
    --             go (Just res) (v:vs) = do
    --                 printf "d_add(%s, %s, %s);\n" res res v :: IO ()
    --                 go (Just res) vs

    --         result <- go Nothing terms
    --         printf "double %s = %s;\n" resultVar result

    --         putStrLn "// ====="

    return ()

{-
main :: IO ()
main = do
    ls <- lines <$> getContents
    ls <- return $ filter (not . null) ls
    let [fs, us, es] = splitOn ["-----"] ls

    -- let funcs = map parseFunc fs
    -- forM_ funcs print

    let funcs = map (either (\e -> error $ "parse failed: " ++ show e) id . parse parseFunction "") fs
    -- forM_ funcs print
    let funcDict = [(name, f) | f@(Function name _ _) <- funcs]

    let bodies = map (fromRight (error "parse failed") . parse parseAssign "") es

    -- mapM_ (putStrLn . ppStmt) bodies

    putStrLn "Constant Propergation"
    bodies <- return $ map (mapStmt constProp) bodies

    putStrLn "Weaken"
    bodies <- return $ map (mapStmt weaken) bodies
    -- bodies <- return $ map (mapStmt removeDiv) bodies

    putStrLn "Inlining"
    bodies <- return $ map (mapStmt $ inlining funcDict) bodies

    putStrLn "Constant Propergation"
    bodies <- return $ map (mapStmt constProp) bodies

    putStrLn "Weaken"
    bodies <- return $ map (mapStmt weaken) bodies
    -- bodies <- return $ map (mapStmt removeDiv) bodies

    putStrLn "Constant Propergation"
    bodies <- return $ map (mapStmt constProp) bodies

    -- mapM_ (putStrLn . ppStmt) bodies

    -- putStrLn "Counting Subexpressions..."
    -- let subexprs = Map.toList $ foldl mergeSubs Map.empty $ map (subExprs . getExpr) bodies
    -- subexprs <- return $ sortOn snd subexprs

    -- forM_ subexprs $ \(e, cnt) -> do
    --     let isTribial =
    --             case e of
    --                 Const _ -> True
    --                 VarExpr (VarRef _) -> True
    --                 _ -> False
    --     when (cnt > 1 && not isTribial) $ do
    --         printf "%d: %s\n" cnt $ ppExpr e

    putStrLn "Lowering"
    let insts = doLower bodies

    forM_ insts $ \(id, inst) -> do
        if id >= 0
            then printf "let e%d = %s\n" id (show inst)
            else printf "%s\n" (show inst)

    putStrLn "Inspecting variable usage..."

    let f :: Set.Set Int -> [(Int, LInst)] -> IO [Int]
        f ss [] = do
            assert (Set.size ss == 0) $ return []

        f ss ((ix, inst): insts) = do
            let refs = instRefs inst
            let nss = Set.delete ix $ foldl (\ss e -> Set.insert e ss) ss refs
            let size = Set.size nss
            -- printf "Live vars: %d\n" size
            rest <- f nss insts
            return $ size : rest

    lives <- f Set.empty (reverse insts)
    printf "Maximum live vars: %d\n" $ maximum lives

    let irFileName = "ir.json"
    printf "Writing IR to %s...\n" irFileName

    let ops = map (\(ix, inst) ->ppInst ix inst) insts
    writeFile irFileName $ "[\n" ++ intercalate ",\n" ops ++ "\n]\n"

    -- putStrLn "Inspecting expressions..."

    -- let f :: [String] -> [(Maybe String, [String])] -> IO [Int]
    --     f _ [] = return []
    --     f vars ((var, _):es) = do
    --         let newVars = maybeToList var ++ vars
    --             restRefs = concat $ map snd es
    --             restVars = filter (`elem` restRefs) newVars
    --             len = length restVars
    --         printf "Live vars: %d\n" len
    --         rest <- f restVars es
    --         return $ len : rest

    -- lives <- f [] bodies

    -- printf "Maximum live vars: %d\n" $ maximum lives

    insts <- codeGen insts

    return ()
-}
