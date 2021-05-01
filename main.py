"""
Ron Mansharof  208839787
Benny Shalom 203500780
"""
def print_Matrix(M):
    '''

    :param M:Matrix
    :return: print matrix
    '''
    for i in range(len(M)):
        print(M[i])
#/////////////////////////////////////////////////
def add_matrix(A,B):
    '''

    :param A:matrix
    :param B: matrix
    :return: add of 2 maricses (A+B)
    '''
    result=copyMatrix(A)
    for i in range(len(A)):
        for j in range(len(A[0])):
            result[i][j] = A[i][j] + B[i][j]
    return result
#/////////////////////////////////////////////////
def zeros_matrix(A):
    '''

    :param A:matrix
    :return: zero matrix in A size
    '''
    b=copyMatrix(A)
    for i in range(len(A)):
        for j in range(len(A[0])):
            b[i][j]=0.0
    return b
#/////////////////////////////////////////////////
def mul_Matrix(A,B):
    '''

    :param A:matrix
    :param B: matrix
    :return: mul A and B (A is in the left, B is in the right)
    '''
    rowsA = len(A)
    colsA = len(A[0])
    rowsB = len(B)
    colsB = len(B[0])
    if colsA != rowsB:
        raise ArithmeticError(
            'Number of A columns must equal number of B rows.')
    C = zeros_matrix(B)
    for i in range(rowsA):
        for j in range(colsB):
            total = 0
            for ii in range(colsA):
                total += A[i][ii] * B[ii][j]
            C[i][j] = total

    return C
#/////////////////////////////////////////////////
def mul_Num_Matrix(n,A):
    '''

    :param n:number
    :param A: matrix
    :return: mul n and A
    '''
    size = len(A)
    M=copyMatrix(A)
    for i in range(size):
        for j in range(size):
            if M[i][j]!=0:
                M[i][j]=n*M[i][j]
    return M
#/////////////////////////////////////////////////
def swapRows(A,row1,row2):
    '''

    :param A: matrix
    :param row1: row
    :param row2: row
    :return: swap row1 and row2 and return the new matrix(A is update)
    '''
    temp=A[row1]
    A[row1]=A[row2]
    A[row2]=temp
    return A
#/////////////////////////////////////////////////
def Pivot(A,row):
    '''

    :param A: matrix
    :param row: row index
    :return: make sure all the numbers on the diagonal are the largest in the column(A is update)
    '''
    maximum=abs(A[row][row])
    help=row
    if row!=len(A)-1:
        for i in range(row+1,len(A)):
            if abs(A[i][row])>=maximum:
                maximum=abs(A[i][row])
                help=i
    if help!=row:
        swapRows(A,row,help)
    return A
#/////////////////////////////////////////////////
def copyMatrix(M):
    '''
    :param M: Matrix
    :return: copy of matrix M
    '''
    size=len(M)
    size1=len((M[0]))
    m=[[0 for i in range(size1)] for j in range(size)]
    for i in range(0,size):
        for j in range(0,size1):
                m[i][j]=M[i][j]
    return m

#/////////////////////////////////////////////////
def identity_Matrix(A):
    '''
    :param A:matrix
    :return: identity matrix in size of A matrix
    '''
    size = len(A)
    b= [[0 for i in range(size)] for j in range(size)]
    for i in range(0,size):
        for j in range(0,size):
            if i==j:
                b[j][i]=1
    return b
#/////////////////////////////////////////////////
def inverseMatrix(A):
    '''
    :param A: matrix
    :return: A^-1
    '''
    A2=copyMatrix(A)
    I=identity_Matrix(A)
    for a in range(len(A)):
        div1 = 1.0 / A2[a][a]
        for j in range(len(A)):
            A2[a][j] *= div1
            I[a][j] *= div1
        for i in list(range(len(A)))[0:a] + list(range(len(A)))[a + 1:]:
            div2 = A2[i][a]
            for j in range(len(A)):
                A2[i][j] = A2[i][j] - div2 * A2[a][j]
                I[i][j] = I[i][j] - div2 * I[a][j]
    return I

# /////////////////////////////////////////////////
def cond(A):
    '''

    :param A:matrix
    :return: cond of matrix A
    '''
    size=len(A)
    a=0
    for i in range(size):
        temp=0
        for j in range (size):
            temp=temp+abs(A[i][j])
        if temp>a:
            a=temp
    return a
# /////////////////////////////////////////////////
def jaacobian_Converge(A):
    '''

    :param A: matrix
    :return:  check if ||G|| < 1 (G = (-D^-1)*(L+U) for jaacobian
    '''
    L,D,U=LDU(A)
    G=add_matrix(L,U)
    D=mul_Num_Matrix(-1,inverseMatrix(D))
    G=mul_Matrix(D,G)
    a=cond(G)
    print("||G||= {0}".format(a))
    if a < 1:
        print("Converge")
        return True
    else:
        print("Does not converge")
        return False
# /////////////////////////////////////////////////
def jaacobian_Calculate(A,b):
    '''

    :param A: matrix in any size not only 3X3
    :param b: matrix of solution
    :return: find x (Ax=b) return x calculate in jaacobian calculation
    '''
    if (jaacobian_Converge(A)== False):
        print("The matrix does not converge")
    else:
        x = copyMatrix(b)
        x = zeros_matrix(b)
        x1 = zeros_matrix(x)
        flag = True
        epsilon= 0.00001
        counter = 0
        print("Count    ", end=" ")
        for i in range(len(x)):
            print("var{0}         ".format(i + 1), end=" ")
        print(" ")
        while flag:
            x=copyMatrix(x1)
            p = 0
            print(counter, end=" ")
            while p < len(x):
                print("     ", end=" ")
                print("%0.5f" % x[p][0], end=" ")
                # print("      ",end = " ")
                p += 1
            for i in range(len(A)):
                temp=b[i][0]
                for j in range(len(A)):
                    if i != j:
                        temp=temp-A[i][j]*x[j][0]
                temp=temp/A[i][i]
                x1[i][0]=temp
            flag=abs(x1[0][0]-x[0][0])>epsilon
            counter+=1
            print(" ")
            #print current x,y,z
        print("Solution={", end=" ")
        for i in range(len(x)):
            print("var{0}=".format(i + 1), end=" ")
            print("%0.6f," % x1[i][0], end=" ")
        print("}")
        print(" ")
# /////////////////////////////////////////////////
def seidel_Converge(A):
    '''

       :param A: matrix
       :return:  check if ||G|| < 1 (G = (-(L+D)^-1)*U) for seidel
       '''
    L,D,U=LDU(A)
    G=inverseMatrix(add_matrix(L,D))
    G=mul_Num_Matrix(-1,G)
    G=mul_Matrix(G,U)
    a=cond(G)
    print("||G||= {0}".format(a))
    if a < 1:
        print("Converge")
        return True
    else:
        print("Does not converge")
        return False
# /////////////////////////////////////////////////
def seidel_Calculation(A,b):
    '''

        :param A: matrix in any size not only 3X3
        :param b: matrix of solution
        :return: find x (Ax=b) return x calculate in seidel calculation
        '''
    if (seidel_Converge(A) == False):
        print("The matrix does not converge")
    else:
        x = copyMatrix(b)
        x = zeros_matrix(b)
        x1 = zeros_matrix(x)
        flag = True
        epsilon = 0.00001
        counter = 0
        #print("Count      x            y            z")
        print("Count    ",end=" ")
        for i in range(len(x)):
            print("var{0}         ".format(i+1),end=" ")
        print(" ")
        while flag:
            x = copyMatrix(x1)
            p=0
            print(counter,end = " ")
            while p<len(x):
                print("     ", end=" ")
                print("%0.6f"%x[p][0],end = " ")
                #print("      ",end = " ")
                p+=1
            for i in range(len(A)):
                temp = b[i][0]
                for j in range(len(A)):
                    if i != j:
                        temp = temp - A[i][j] * x1[j][0]
                temp = temp / A[i][i]
                x1[i][0] = temp
            flag = abs(x1[0][0] - x[0][0]) > epsilon
            counter += 1
            print(" ")
        print("Solution={", end=" ")
        for i in range(len(x)):
            print("var{0}=".format(i+1), end=" ")
            print("%0.6f,"%x1[i][0],end = " ")
        print("}")
        print(" ")

# /////////////////////////////////////////////////
def LDU(A):
    '''

    :param A:matrix
    :return: L , D , U matrices
    '''
    L=copyMatrix(A)
    D=copyMatrix(A)
    U=copyMatrix(A)
    size=len(A)
    for i in range(size):
        for j in range(size):
            if i != j:
                D[i][j]=0
    for i in range(size):
        L[i][i]=0
        U[i][i]=0
    for i in range(1,size):
        for j in range(i-1,-1,-1):
            L[j][i]=0
    for i in range(size):
        for j in range(i,size):
            U[j][i]=0
    return L,D,U
# /////////////////////////////////////////////////
def menu(A,b):
    t=int(input("Press 1 for jaacobian / 2 for seidel / 3 for both\n"))

    print(t)
    if t == 1:
        print("Jaacobian:")
        jaacobian_Calculate(A,b)
    elif t==2:
        print("Seidel:")
        seidel_Calculation(A,b)
    elif t == 3:
        print("Jaacobian:")
        jaacobian_Calculate(A, b)
        print("Seidel:")
        seidel_Calculation(A, b)
    else:
        print("Eror")
# /////////////////////////////////////////////////
A=[[4,2,0],[2,10,4],[0,4,5]]
b=[[2],[6],[5]]
menu(A,b)
