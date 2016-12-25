sage_server.MAX_STDOUT_SIZE=sage_server.MAX_STDOUT_SIZE*10

#---------------------basic parameters
q=31
v=56
d=2
op=7
o=d*op
n=v+o
K.<a>=GF(q)
P=PolynomialRing(K,'x',n)
x=P.gens()
x_vec=vector(x)



#---------------------generating public key
A00=[0 for i in range(o)]
A01=[0 for i in range(o)]
B0=[0 for i in range(o)]
B=[[] for i in range(o)]
A=[0 for i in range(o)]
D=[0 for i in range(o)]
D[0]=random_matrix(K,v,v)

rotate_matrix=[]
I=matrix.identity(K,v,v)
rotate_matrix.append(I.row(v-1).list())
for i in range(v-1):
    rotate_matrix.append(I.row(i).list())
rotate_matrix=matrix(K,v,v,rotate_matrix)



for i in range(o):#-------------------B[i] is vv part
    D[i]=rotate_matrix^i*D[0]

Ar=[0 for i in range(op)]
for r in range(op):
    A01[r]=[]
    for i in range(v*op):
        A01[r].append(K.random_element())
    Ar[r]=matrix(K,op,v,A01[r]).transpose()
    for i in range(v*(o-op)):
        A01[r].append(K(0))
    A01[r]=matrix(K,o,v,A01[r]).transpose()

rotate_matrix=[]
for i in range(op):
    for j in range(o):
        if(j==(o-op+i)):
            rotate_matrix.append(K(1))
        else:
            rotate_matrix.append(K(0))
for i in range(o-op):
    for j in range(o):
        if(j==i):
            rotate_matrix.append(K(1))
        else:
            rotate_matrix.append(K(0))
rotate_matrix=matrix(K,o,o,rotate_matrix).transpose()
for i in range(o):
    A01[i]=A01[i%op]*rotate_matrix^(i//op)

for i in range(o):
    A00[i]=D[i]
for i in range(o):
    A[i]=matrix.block([[A00[i],A01[i]],[matrix(K,o,v),matrix(K,o,o)]])

#---------------omit linear and constant parts for simplicity

F=[0 for i in range(o)]
for i in range(o):
    F[i]=A[i]#-------------------central quadratic matrix
while true:
    T=random_matrix(K,n)
    if T.is_invertible():
        break
Tt=T.transpose()
P=[0 for i in range(o)]
for i in range(o):
    P[i]=Tt*F[i]*T#------------public key of MB_UOV



#-------------constructing equivalent key
o1=op
o2=op
T1=T.submatrix(0,0,v,v);
T2=T.submatrix(0,v,v,o1);
T3=T.submatrix(0,v+o1,v,o2);
T4=T.submatrix(v,0,o1,v);
T5=T.submatrix(v,v,o1,o1);
T6=T.submatrix(v,v+o1,o1,o2);
T7=T.submatrix(v+o1,0,o2,v);
T8=T.submatrix(v+o1,v,o2,o1);
T9=T.submatrix(v+o1,v+o1,o2,o2)
O1=T1
O4=T4
O7=T7
Tp1=T1.inverse()*T2
Tp2=T1.inverse()*T3
O5=T5-T4*T1.inverse()*T2
Tp3=O5.inverse()*(T6-T4*T1.inverse()*T3)
O8=T8-T7*T1.inverse()*T2
O9=T9-T7*T1.inverse()*T3-O8*Tp3
Om=block_matrix(3,3,[[O1,matrix(K,v,o1),matrix(K,v,o2)],[O4,O5,matrix(K,o1,o2)],[O7,O8,O9]])
Tp=block_matrix(3,3,[[identity_matrix(K,v),Tp1,Tp2],[matrix(K,o1,v),identity_matrix(K,o1),Tp3],[matrix(K,o2,v),matrix(K,o2,o1),identity_matrix(K,o2)]])#------Equivalent key
#--------every layer of the Equivalent key
Ti=[0 for i in range(n+1)]
for i in range(v+1,v+op+1,1):
    Ti[i]=identity_matrix(K,n)
    Ti[i].set_block(0,i-1,matrix(n,1,Tp.column(i-1).list()))
Tx=identity_matrix(K,n)
Tx.set_block(0,v,Tp.submatrix(0,v,v,op))
Ty=Tx.inverse()*Tp

for i in range(v+op+1,n+1,1):
    Ti[i]=identity_matrix(K,n)
    Ti[i].set_block(0,i-1,matrix(n,1,Ty.column(i-1).list()))

#------------------------------  F \circ T = F \circ Om \circ Om^(-1) \circ T
#-------------the new central matrix is F \circ Om
#-------------the new affine transformation is Om^(-1)\circ T, it can be reconstructed by Rainbow-Band-Separation attack
#-------------Correctness check
for i in range(o):
    print i,"-th polynomial is equivalent:",(Om*Tp).transpose()*F[i]*(Om*Tp)==P[i]
print "new affine transformation:";
print Tp.str()


for i in range(o):
    print i,"-th new central matrix:"
    print (Om.transpose()*F[i]*Om).str()
    print " "





