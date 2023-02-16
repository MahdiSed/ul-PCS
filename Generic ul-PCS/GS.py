from charm.toolbox.pairinggroup import G1,pair
from charm.toolbox.secretutil import SecretUtil

class GS():
    def __init__(self, groupObj):
        global util, group
        util = SecretUtil(groupObj)        
        group = groupObj

    def MatMul(self,a,b):
            c = []
            for i in range(len(a)):
                s = 1
                for j in range(len(b[0])):
                    s *= a[i]**b[j][i]
                c.append(s)
            return c
    
    def sampleParams(self,pp):
        rho, zeta = group.random(), group.random()
        sigma, omega = group.random(), group.random()
        vv1 = [pp['G1']**zeta, pp['G1']]
        vv2 = [pp['G2']**omega, pp['G2']]
        ww1 = [a**rho for a in vv1]
        ww2 = [a**sigma for a in vv2]
        uu1 = [0, ww1[1]*pp['G1']]
        uu2 = [0, ww2[1]*pp['G2']]
        Zeta = [-zeta**(-1), 1]
        Omega = [-omega**(-1), 1]
        ck = {'uu1':uu1, 'vv1':vv1, 'ww1':ww1, 'uu2':uu2, 'vv2':vv2, 'ww2':ww2}
        #xk = {'ck':ck, 'Zeta':Zeta, 'Omega': Omega}
        return ck


    def ParamGen(self,X,Y,c_a,c_b):
        m = len(c_a); n = len(c_b)
        R_x = [group.random() for _ in range(m)]
        R_y = [group.random() for _ in range(n)]
        S_x = [0 for _ in range(m)]
        S_y = [0 for _ in range(n)]
        return (R_x,R_y,S_x,S_y)
    
    def commit(self, ck ,X, Y, R_x,R_y,S_x,S_y):
        cc_x={}; cc_y={}
        for i in range(len(X)):
            cc_x[i] = [(ck['vv1'][0]**R_x[i])*(ck['ww1'][0]**S_x[i]),\
                        X[i]*(ck['vv1'][1]**R_x[i])*(ck['ww1'][1]**S_x[i])]
        for i in range(len(Y)):
            cc_y[i] = [(ck['vv2'][0]**R_y[i])*(ck['ww2'][0]**S_y[i]),\
                        Y[i]*(ck['vv2'][1]**R_y[i])*(ck['ww2'][1]**S_y[i])]
        return cc_x, cc_y 

    def prove(self,ck, X, Y, R_x,R_y,S_x,S_y, Gamma, cc_x, cc_y):
        alpha, beta, gamma, delta = group.random(), group.random(),group.random(),group.random()
        pi_v1={}; pi_w1={}; pi_v2={}; pi_w2={}
        for i in range(len(X)):
            aux_R = GS.MatMul(self,R_x, [*zip(*Gamma)])
            #aux_S = GS.MatMul(self,S_x, [*zip(*Gamma)])
            for j in range(2):
                pi_v1[i,j] = (cc_y[i][j]**aux_R[j])*(ck['vv2'][j]**alpha)*(ck['ww2'][j]**beta)
                pi_w1[i,j] = (cc_y[i][j]**S_x[j])*(ck['vv2'][j]**gamma)*(ck['ww2'][j]**delta)
            aux_R = GS.MatMul(self,R_y, [*zip(*Gamma)])
            #aux_S = GS.MatMul(self,S_y, [*zip(*Gamma)])
            pi_v2[i,0] = group.init(0,G1); pi_v2[i,1] = (X[i]**aux_R[1])/((ck['vv1'][1]**alpha)*(ck['ww1'][1]**gamma))
            pi_w2[i,0] = group.init(0,G1); pi_w2[i,1] = (X[i]**S_y[1])/((ck['vv1'][1]**beta)*(ck['ww1'][1]**delta))
        pi = {'pi_v1':pi_v1, 'pi_w1':pi_w1, 'pi_v2':pi_v2, "pi_w2":pi_w2}
        return pi
    
    def verifyProof(self,pp,ck,Pi,CC_x,CC_y,Gamma):
        LHS=1; RHS=1; p1={}; p2={}; aux={}
        
        result=True
        for ctr in range(1,len(Gamma)+1):
            cc_x = CC_x[ctr]; cc_y = CC_y[ctr]; pi=Pi[ctr]; gamma=Gamma[ctr]
            m=len(cc_x); n=len(cc_y)
            for vv1 in range(2):
                for vv2 in range(2):
                    for i in range(m):
                        p1[i] = cc_x[i][vv1]
                        p2[i] = 1
                        for j in range(n):
                            p2[i] *= cc_y[j][vv2]**gamma[j][i]
                    for i in range(2):
                            p1[m+i] = ck['vv1'][vv1]**(-1)
                            p2[m+i] = pi['pi_v2'][i,vv2]
                            p1[m+2+i] = pi['pi_v1'][i,vv1]
                            p2[m+2+i] = ck['vv2'][vv2]**(-1)
                    for i in range(m):
                        LHS *= pair(p1[i],p2[i])
                    if LHS==1:
                        result=False
        return result
