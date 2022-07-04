#https://www.nature.com/articles/s41598-021-85434-9.pdf
import numpy as np

#Initial data
r_1 = #input radius
x_p = # input position
Q = #input blood flow
Omega = #Vascular domain
N_TF = #Number of terminals in the resulting tree
f_r = #correction step factor (0,1)

def optimazation():
    #Define triangle (x^d_i, x^p_j, x^d_j)
    dv = #discretazion between x_j^p and x_j^d


lc = np.sqrt(Omega_size/np.pi) #equation (13) D = 2
lmin = lc*np.sqrt(nu/(Nt+1)) #equation (12) nu = tuning parameter , Nt = number of terminal vessels
#note, there is an update rule/failsafe on lmin, should be considered later


while l[i]/r[i] <= 2 and l[i] <= l_min:
    x_d = #randomly generate point inside omega (uniform distrubution)
    v[0] = (r_1, x_p, x_d)

# Add v[0] to F (F = current tree)
for i in range(1, N_TF):
    L = []
    while np.size(L) == 0:
        while np.abs(x_d[i]-x_n) <= lmin:
            if iter == N_fail:
                lmin = f_r*lmin
            x_d = #randomly generate point inside omega (uniform distrubution)
            #determine x_n as the nearest point of F to x_d[i]
        #Determine the subset Fn of possible vessels to connect x_d[i]
        for v[j] in Fn:


import numpy as np

l=5000;
a=zeros(l,2);
for i=2:l
    b=200*(rand(1,2)-.5);
    [~,w]=min(sum((a-b).^2,2));
    n=a(w,:);
    k=atan2(b(2)-n(2),b(1)-n(1));
    c=[n(1)+cos(k)  n(2)+sin(k)];
    a(i,:)=c;
    plot([n(1) c(1)],[n(2) c(2)],'g','Linew',2)
    hold on
end
