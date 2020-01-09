
import matplotlib.pyplot as plt

plt.hist(mye,bins=50,alpha=0.75,normed=True,label='Prob having e > 0 for circular orbit')
plt.hist(e_vec,bins=50,alpha=0.75,normed=True,label='PDF for e coming from MCMC')
plt.legend(loc=0, ncol=1,scatterpoints=1,numpoints=1,frameon=False)
plt.xlabel('eccentricity')
plt.ylabel('Probability')
plt.savefig('probs.png')
plt.show()
