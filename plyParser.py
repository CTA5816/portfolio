##adapted from ply example calc.py

import expression as xp
import operation as op
import re
import expToSurface
import pygame

def mkstr_plus(opList):
    opStr = opList[0].getString() + "+" + opList[1].getString()
    return opStr

def mkstr_minus(opList):
    opStr = opList[0].getString() + "-" + opList[1].getString()
    return opStr
	
def mkstr_times(opList):
    opStr = opList[0].getString() + '\u00B7' + opList[1].getString()
    return opStr

def mkstr_silenttimes(opList):
    opStr = opList[0].getString() + opList[1].getString()
    return opStr

def mkstr_slashdiv(opList):
    opStr = opList[0].getString() + "/" + opList[1].getString()
    return opStr
	
def mkstr_caretpow(opList):
    opStr = opList[0].getString() + "^" + opList[1].getString()
    return opStr

def mkstr_fullprns(opList):
    opStr = "(" + opList[0].getString() + ")"
    return opStr

def mkstr_lprns(opList):
    opStr = "(" + opList[0].getString()
    return opStr

def mkstr_rprns(opList):
    opStr = opList[0].getString() + ")"
    return opStr

def mkstr_silentprns(opList):
    opStr = opList[0].getString()
    return opStr

def mkstr_error(opList):
    return "ERROR"
	
ADD_OP = op.Operator('+', mkstr_plus)
SUB_OP = op.Operator('-', mkstr_minus)
MUL_OP = op.Operator('\u00B7', mkstr_times)
SMUL_OP = op.Operator('*', mkstr_silenttimes)
DIV_OP = op.Operator('/', mkstr_slashdiv)
POW_OP = op.Operator('^', mkstr_caretpow)
FPRN_OP = op.Operator('()', mkstr_fullprns)
LPRN_OP = op.Operator('(', mkstr_lprns)
RPRN_OP = op.Operator(')', mkstr_rprns)
SPRN_OP = op.Operator('{}', mkstr_silentprns)

ERR_OP = op.Operator("ERROR", mkstr_error)

tokens = (
    'VAR', 'NUMBER',
    '1L1R_OP_L0', '1L1R_OP_L1', '1L1R_OP_R2',
    'LPAREN','RPAREN', 'ULPAREN', 'URPAREN', 'LBRACK', 'RBRACK', 'NUM_CURSOR', 'VAR_CURSOR', 'UNF_CURSOR', 'LPRN_CURSOR', 'RPRN_CURSOR', 'ULPRN_CURSOR', 'URPRN_CURSOR', 'LBRK_CURSOR', 'RBRK_CURSOR'
    )

# Tokens
#t_EQUALS  = r'='
t_LPAREN  = r'\('
t_RPAREN  = r'\)'
t_ULPAREN = r'\u2985'
t_URPAREN = r'\u2986'
t_LBRACK  = r'\{'
t_RBRACK  = r'\}'
#t_NAME    = r'[a-zA-Z_][a-zA-Z0-9_]*'
t_UNF_CURSOR  = r'((?<=[^a-zA-Z\d\)\u2986\}])|(?<=^))\|(?=[^a-zA-Z\d\(\u2985\{]|$)'
#t_PRN_CURSOR = r'\|(?=\()|(?<=\))\|'
t_LPRN_CURSOR = r'\|(?=\()'
t_RPRN_CURSOR = r'(?<=\))\|'
t_ULPRN_CURSOR = r'\|(?=\u2985)'
t_URPRN_CURSOR = r'(?<=\u2986)\|'
t_LBRK_CURSOR = r'\|(?=\{)'
t_RBRK_CURSOR = r'(?<=\})\|'

def t_1L1R_OP_L0(t):
    r'\+|-'
    if t.value == "+":
            t.value = ADD_OP
    elif t.value == '-':
            t.value = SUB_OP
    else:
            t.value = ERR_OP
    return t

def t_1L1R_OP_L1(t):
    r'\*|/|\u00B7'
    if t.value == '\u00B7':
            t.value = MUL_OP
    elif t.value == '*':
            t.value = SMUL_OP
    elif t.value == '/':
            t.value = DIV_OP
    else:
            t.value = ERR_OP
    return t	

def t_1L1R_OP_R2(t):
    r'\^'
    if t.value == '^':
            t.value = POW_OP
    else:
            t.value = ERR_OP
    return t
    
def t_NUM_CURSOR(t):
    r'(\d+\|\d*)|(\d*\|\d+)'
    numString = t.value
    cursor_idx = numString.index('|')
    numString = numString.replace('|', '')
    t.value = xp.NoOpExpression(numString, True, cursor_idx)
    return t
    
def t_VAR_CURSOR(t):
    r'([a-zA-Z]\|)|(\|[a-zA-Z])'
    varString = t.value
    cursor_idx = varString.index('|')
    varString = varString.replace('|', '')
    t.value = xp.NoOpExpression(varString, True, cursor_idx)
    return t
	
def t_VAR(t):
    r'[a-zA-Z]'
    t.value = xp.NoOpExpression(t.value)
    return t
	
def t_NUMBER(t):
    r'\d+(\.\d+)?'
    t.value = xp.NoOpExpression(str(t.value))
    return t

# Ignored characters
t_ignore = " "

def t_newline(t):
    r'\n+'
    t.lexer.lineno += t.value.count("\n")
    
def t_error(t):
    global error
    error=True
    print("Illegal character '%s'" % t.value[0])
    t.lexer.skip(1)
    
# Build the lexer
import ply.lex as lex
lexer = lex.lex(reflags=re.UNICODE)

# Parsing rules

#def p_statement_assign(p):
#    'statement : statement EQUALS exp0'
#    

def p_statement_exp0(p):
    'statement : exp0'
    p[0] = p[1]

def p_exp0_exp1(p):
    'exp0 : exp1'
    p[0] = p[1]

def p_exp1_exp2(p):
    'exp1 : exp2'
    p[0] = p[1]	

def p_exp2_exp3(p):
    'exp2 : exp3'
    p[0] = p[1]	

def p_exp3_parens(p):
    'exp3 : LPAREN exp0 RPAREN'
    p[0] = xp.Expression(FPRN_OP, [p[2]])

def p_exp3_clparens(p):
    'exp3 : LPRN_CURSOR LPAREN exp0 RPAREN'
    p[0] = xp.Expression(FPRN_OP, [p[3]]).addCursor(0)

def p_exp3_crparens(p):
    'exp3 : LPAREN exp0 RPAREN RPRN_CURSOR'
    p[0] = xp.Expression(FPRN_OP, [p[2]]).addCursor(1)

def p_exp3_ulparens(p):
    'exp3 : LPAREN exp0 URPAREN'
    p[0] = xp.Expression(LPRN_OP, [p[2]])

def p_exp3_culparens(p):
    'exp3 : LPRN_CURSOR LPAREN exp0 URPAREN'
    p[0] = xp.Expression(LPRN_OP, [p[3]]).addCursor(0)

def p_exp3_cluparens(p):
    'exp3 : ULPRN_CURSOR ULPAREN exp0 RPAREN'
    p[0] = xp.Expression(RPRN_OP, [p[3]]).addCursor(0)

def p_exp3_urparens(p):
    'exp3 : ULPAREN exp0 RPAREN'
    p[0] = xp.Expression(RPRN_OP, [p[2]])

def p_exp3_curparens(p):
    'exp3 : ULPAREN exp0 RPAREN RPRN_CURSOR'
    p[0] = xp.Expression(RPRN_OP, [p[2]]).addCursor(1)	

def p_exp3_cruparens(p):
    'exp3 : LPAREN exp0 URPAREN URPRN_CURSOR'
    p[0] = xp.Expression(LPRN_OP, [p[2]]).addCursor(1)	

def p_exp3_bracks(p):
    'exp3 : LBRACK exp0 RBRACK'
    p[0] = xp.Expression(SPRN_OP, [p[2]])

def p_exp3_clbracks(p):
    'exp3 : LBRK_CURSOR LBRACK exp0 RBRACK'
    p[0] = xp.Expression(SPRN_OP, [p[3]]).addCursor(0)

def p_exp3_crbracks(p):
    'exp3 : LBRACK exp0 RBRACK RBRK_CURSOR'
    p[0] = xp.Expression(SPRN_OP, [p[2]]).addCursor(1)

def p_exp3_number(p):
    'exp3 : NUMBER'
    p[0] = p[1]

def p_exp3_cnumber(p):
    'exp3 : NUM_CURSOR'
    p[0] = p[1]

def p_exp3_cvar(p):
    'exp3 : VAR_CURSOR'
    p[0] = p[1]

def p_exp3_var(p):
    'exp3 : VAR'
    p[0] = p[1]

def p_exp3_empty(p):
    'exp3 : empty'
    p[0] = xp.NoOpExpression("")

def p_exp3_cempty(p):
	'exp3 : UNF_CURSOR empty'
	p[0] = xp.NoOpExpression("", True, 0)

def p_empty(p):
     'empty :'
     pass

def p_exp0_exp0ops(p):
    'exp0 : exp0 1L1R_OP_L0 exp1'
    expList = [p[1], p[3]]
    p[0] = xp.Expression(p[2], expList)

def p_exp1_exp1ops(p):
    'exp1 : exp1 1L1R_OP_L1 exp2'
#    if p[2].strrep == '/':
#		
#    else:
#        expList = [p[1], p[3]]
    expList = [p[1], p[3]]
    p[0] = xp.Expression(p[2], expList)

def p_exp1_exp2ops(p):
    'exp2 : exp3 1L1R_OP_R2 exp2'
    expList = [p[1], p[3]]
    p[0] = xp.Expression(p[2], expList)

def p_error(t):
    global error
    error=True
    print("Syntax error at '%s'" % t.value)

import ply.yacc as yacc
parser = yacc.yacc()

def preprocess(inputStr):
    output_str = inputStr	
    #Insert implicit multiplication
    output_str = re.sub('(?<=\w|\))(?=\|?\()|(?<=\))(?=\|?\w)', '*', output_str)
    output_str = re.sub('(?<=\w|\})(?=\|?\{)|(?<=\})(?=\|?\w)|(?<=\d|[a-zA-Z])(?=\|?[a-zA-Z])|(?<=[a-zA-Z])(?=\|?\d)', '*', output_str)
    return output_str

error=False
def get_exp(inputStr):
    if inputStr[1:]=='loading':
        return xp.Expression(op.Operator('loading...',lambda x:''),[xp.NoOpExpression("")])
    global error
    error=False
    print('processing',inputStr.replace('\u2986','$').replace('\u2985','#'))
    inputStr = preprocess(inputStr)
    print('got',inputStr.replace('\u2986','$').replace('\u2985','#'))
    resultingExpression = parser.parse(inputStr)
    if error:
        raise ValueError("Expression could not parse correctly.")
    else:
        return resultingExpression
        #string = resultingExpression.getString()
        #Surface = expToSurface.smartSurface(string)
        #return Surface.surface