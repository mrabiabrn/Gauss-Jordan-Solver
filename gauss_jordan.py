
class GaussJordanSolver(object):
    n = ""
    matrix = ""
    original = ""

    def __init__(self, n, matrix):
        self.n = n
        self.matrix = matrix
        self.original = matrix[:]

    def multiply_row(self, row, cons):
        self.matrix[row] = [ i/cons if abs(i/cons) > 1e-7 else 0 for i in self.matrix[row]]

    def interchange(self, i1, i2):
        self.matrix[i1], self.matrix[i2] = self.matrix[i2], self.matrix[i1]

    @staticmethod 
    def get_determinant(grid):
        if len(grid) == 2:
            return grid[0][0]*grid[1][1]-grid[0][1]*grid[1][0]
        det = 0
        mult = 1
        for i in range(len(grid)):
            num = grid[0][i]
            grid2 = [row[0:i] + row[i+1:len(grid)] for row in grid[1:] ]
            det2 = GaussJordanSolver.get_determinant(grid2)
            det = det + num * mult * det2
            mult = mult * -1 
        
        return det

    def inverse(self):
        res = []
        det = 0

        A = [[row[i] for i in range(len(row)) if i != self.n ] for row in self.original]

        for i in range(len(A)):
            temp = []
            for j in range(len(A)):
                matrix2 = [row[0:j] + row[j+1:len(A)] for row in A[0:i] + A[i+1:len(A)]]
                det2 = GaussJordanSolver.get_determinant(matrix2)  
                temp.append(det2 * (-1)**(i+j))
                if(i==0):
                    det = det + A[i][j] * det2 * (-1)**j
                
            res.append(temp)
        
        res = list(map(list,zip(*res)))

        for i in range(len(res)):
            for j in range(len(res)):
                res[i][j] = res[i][j]/det + 0
        
        return res


    # matrix[i1] + k*matrix[i2] -> matrix[i1]
    def add(self, i1, i2, k): 
        row1, row2 = self.matrix[i1], self.matrix[i2]  
        self.matrix[i1] = [r1 + k*r2 if abs(r1 + k*r2) > 1e-7 else 0 for (r1, r2) in zip(row1, row2)]

    # returns the id of the next row, -1 if all zero
    def next_row(self, rowStart, columnIndx):

        for i in range(rowStart, self.n):
            if self.matrix[i][columnIndx] != 0:
                return i
            
        return -1

    def find_rank(self):
        i,row,rank,last = 0, 0, 0, 0

        while(i<=self.n and row<self.n):
            if(self.matrix[row][i] == 1 ):
                rank = rank + 1
                last = i
                i = i + 1
                row = row + 1
            elif(i==self.n):
                i = last + 1
                row = row + 1
            else:
                i = i + 1

        return rank-1 if last==self.n else rank, rank

    @staticmethod
    def print_matrix(matrix):
        stri = ""
        for i in range(len(matrix)):
            if( i!=0):
                stri = stri + " "*12
            for j in range(len(matrix[0])):
                num = str(round(matrix[i][j],3))
                stri = stri + " "*(7-len(num)) + num
            stri = stri + "\n"
 
        return stri
            

    # rank1 = rank(A), rank2 = rank(A|b),
    def find_solution(self, rank1, rank2):

        solution = []

        # unique solution
        if rank1 == self.n:

            for row in self.matrix:
                solution.append(row[-1])
            
            
            return  "Unique solution: " + str(solution)[1:-1] + "\n" + "Inverted A: " + GaussJordanSolver.print_matrix(self.inverse())

        # infinetely many solution
        elif rank1 == rank2:

            variables = []  # non-basic variables (arbitrary)
            count = 0       # num of variables whose values are found    
            for indx in range(rank1):
                cur_row = self.matrix[indx]
                for i in range(count,self.n ):
                    cur_num = cur_row[i] 
                    # non-basic
                    if cur_num == 0:
                        solution.append(0)
                        count +=1
                        variables.append("x" + str(count))
                    # basic (pivot)
                    elif cur_num == 1:
                        solution.append(cur_row[-1])
                        count+=1
                        break
            
            solution.extend([0]*(self.n-count))
            variables.extend(["x" + str(i) for i in range(count,n)])

            return  str(["x" + str(i) for i in range(n)])[1:-1] + " are variables." + "\n" + "Arbitrary variables: " + str(variables)[1:-1] + "\n" +"Arbitrary solution: " + str(solution)[1:-1]
            
        # no solution
        else:
            return "Inconsistent problem"

    def solver(self):

        rowStart = 0
        columnIndx = 0

        while(rowStart<self.n and columnIndx<self.n+1):  #
            row = self.next_row(rowStart, columnIndx)
           
            if(row == -1):
                columnIndx = columnIndx + 1
                continue

            self.interchange(rowStart,row)

            num = self.matrix[rowStart][columnIndx]
            self.multiply_row(rowStart, num)

            for j in range(self.n):
                if(j == rowStart):
                    continue
                else:
                    num2 = self.matrix[j][columnIndx]
                    self.add(j, rowStart, -num2)
            
            columnIndx = columnIndx + 1
            rowStart = rowStart + 1

        
        rank1, rank2 = self.find_rank()
        print(self.find_solution(rank1,rank2))


def read_file(file_name):
    f = open(file_name, "r")
    data = f.read().split()
    matrix = []
    n = int(data[0])
    for i in range(n):
        temp = []
        for j in range(n+1):
            temp.append(float(data[i*(n+1) + j + 1]))
        matrix.append(temp)
    f.close()
    return n , matrix



n, matrix = read_file("Data.txt")
g1 = GaussJordanSolver(n, matrix)
g1.solver()
    
