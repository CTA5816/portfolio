import operation as opn
import numpy as np
import wolframalpha

APP_ID = "8T8YA5-3V337TXULH"
client = wolframalpha.Client(APP_ID)

class Expression:
    #Initializer for composite expression
    def __init__(self, op, expList, cursor=False, cursor_idx=None):
        '''
        Initializes a composite expression
        Parameters:
            self: Object being initialized
            op: operation the expression contains
            explist: List of expressions the expression contains (assumed to be ordered properly; for example, it would be [2,1] for 2+1 and [1,2] for 1+2)
            parens: Whether or not this expression has parentheses
            parens_solid: Encodes graphical information about the solidity (white-black ratio) of parentheses if the expression has them. [1,0] means the expression's left paren is black and right paren is grey, for example
        '''
        self.expList = expList #Stores expressions that this expression comprises
        self.op = op
        self.cursor = cursor
        self.cursor_idx = cursor_idx
        self.parent = None
    #Gets string representation of expression
    def getString(self):
        if type(self)==NoOpExpression:
            return self.strRep
        returnString = self.op.makeStr(self.expList)
        return returnString

    #Evaluates expression
    def eval(self):
        expString = self.getString()
        res = client.query(expString) #Gets result from WolframAlpha
        return next(res.results).text

    #Adds cursor and index
    def addCursor(self, idx):
        self.cursor = True
        self.cursor_idx = idx
        return self

    def __repr__(self):
        return repr(self.op)+' of ('+') and ('.join([repr(k) for k in self.expList])+')'

    def assign_parents(self):
        currentNode = self
        if type(currentNode) == Expression:
            for node in currentNode.expList:
                node.parent = currentNode
            for node in currentNode.expList:
                node.assign_parents()




class NoOpExpression(Expression):
    #Initializer for a no-operator expression
    def __init__(self, strRep, cursor=False, cursor_idx=None):
        '''
        Parameters:
            self: Object being initialized
            strRep: String representation of expression.
            '''
        self.op = opn.Operator('',lambda x:''.join(x)) #null operator
        self.strRep = strRep
        self.cursor = cursor
        self.cursor_idx = cursor_idx
        self.parent = None

	#Adds cursor and index
    def addCursor(self, idx):
        self.cursor = True
        self.cursor_idx = idx
		
    def __repr__(self):
        return self.strRep

    def assign_parents(self):
        return 0


def compare_exp(first,second):
    if not type(first) == type(second):
        return False
    if type(first) == Expression:
        if not first.op == second.op:
            return False
        if not len(first.expList) == len(second.expList):
            return False
        for i,j in zip(first.expList, second.expList):
            if not compare_exp(i,j):
                return False
            return True
    if type(first) == NoOpExpression:
        if not first.strRep == second.strRep:
            return False
        return True

    print("invalid comparison")
    return False

def flatten(exp):
    if type(exp) == NoOpExpression:
        return [exp]
    allexps = []
    allexps += [exp]
    for child in exp.expList:
        allexps += flatten(child)
    return allexps


        

