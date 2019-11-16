#Let us do the plots here
from matplotlib import gridspec
from matplotlib.colors import LogNorm

if ( is_seaborn_plot ):
  import seaborn as sns
  sns.set(style='ticks')
  sns.set_color_codes(seaborn_palette)

fsx = figure_size_x
fsy = figure_size_y
fos = font_size_label

vari = params[0]
posterior = params[1]
chi2 = params[2] + params[3]

#===========================================================
#              plot chains
#===========================================================

def plot_chains():
  plt.xlabel('iteration')
  plt.ylabel('Reduced $\chi^2$')
  plt.hist2d(vari,chi2/(ndata-npars),bins=50)
  fname = outdir+'/'+tic+'_chains.pdf'
  print('Creating ', fname)
  plt.savefig(fname,bbox_inches='tight')
  plt.close()

def plot_postiter():
  plt.xlabel('iteration')
  plt.ylabel('$\ln \mathcal{L}$')
  plt.hist2d(vari,posterior,bins=50)
  fname = outdir+'/'+tic+'_likelihood.pdf'
  print('Creating ', fname)
  plt.savefig(fname,bbox_inches='tight')
  plt.close()

#This is a general routine to plot a two boxes plot
#first plot is a time series with a model and data
#secod plot is a tme series with residuals
def create_nice_plot(mvector,dvector,labels,mlabels,inst_labels,fname,\
                     fsx=fsx,fsy=fsy,plot_residuals=True,
                     model_colors='k',model_alpha=1):

  print('Creating ', fname)
  #
  tmodel  = mvector[0]
  #
  tdata   = dvector[0]
  ydata   = dvector[1]
  edata   = dvector[2]
  ejdata  = dvector[3]
  res     = dvector[4]
  labvec  = dvector[5]

  #
  label_x   = labels[0]
  label_y   = labels[1]
  label_res = labels[2]
  #
  if model_colors.__class__ != list:
    model_colors = [None]*(len(mvector)-1)
    model_colors[:] = str('k')
  if model_alpha.__class__ != list: model_alpha = [1.]*(len(mvector)-1)
  #

  plt.figure(1,figsize=(fsx,fsy))
  if plot_residuals:
    gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[3., 1.])
  else:
    gs = gridspec.GridSpec(nrows=1, ncols=1)
  gs.update(hspace=0.00)
  #Timeseries with data and model plot
  ax0 = plt.subplot(gs[0])
  plt.tick_params(labelsize=fos,direction='in')
  plt.minorticks_on()
  plt.xlabel("")
  plt.ylabel(label_y,fontsize=fos)
  #PLOT MODELS
  for j in range(1,len(mvector)):
    plt.plot(tmodel,mvector[j],linewidth=2.0,label=mlabels[j-1],\
             color=model_colors[j-1],alpha=model_alpha[j-1],zorder=3)
  #PLOT DATA
  for j in range(0,len(tdata)):
    #This one plots the jitter error bars
    plt.errorbar(tdata[j],ydata[j],ejdata[j],fmt='.',\
    alpha=0.9,color='#C0C0C0',markersize=rv_markersize,fillstyle='none',zorder=1)
    #This one plots the legacy error bars
    plt.errorbar(tdata[j],ydata[j],edata[j],\
    fmt=mark[labvec[j]],alpha=1.0,color=rv_colors[labvec[j]], \
    markersize=rv_markersize,fillstyle=rv_fillstyle,zorder=2)
  for j in range(0,len(inst_labels)):
    plt.errorbar(-1,np.mean(ydata),1e-5,fmt=mark[j],label=inst_labels[j],\
    alpha=1.0,color=rv_colors[j], \
    markersize=rv_markersize,zorder=2)
  #
  if ( is_rv_legend ): plt.legend(ncol=1,scatterpoints=1,numpoints=1,frameon=True,fontsize=fos*0.7)
  #
  #plt.xticks(np.arange(0.,1.01,0.1))
  plt.tick_params( axis='x',which='both',direction='in',labelbottom=False)
  plt.tick_params( axis='y',which='both',direction='in')
  if not plot_residuals:
    plt.xlabel(label_x,fontsize=fos)
    plt.tick_params( axis='x',which='both',direction='in',labelbottom=True)
  #
  yylims = ax0.get_ylim()
  miy = int(max(abs(yylims[0]),abs(yylims[1])))
  if ( miy == 0): #To avoid errors for really small RV variations
    miy = max(abs(yylims[0]),abs(yylims[1]))
#  plt.ylim(-miy,miy)
  #
  #if ( select_y_rv ):
  #  plt.ylim(rv_lim_min,rv_lim_max)
  plt.xlim(min(tmodel),max(tmodel))
  #
  #NEW SUBPLOT: RESIDUALS
  #
  if plot_residuals:
    ax1 = plt.subplot(gs[1])
    plt.tick_params(labelsize=fos,direction='in')
    plt.xlabel(label_x,fontsize=fos)
    plt.tick_params( axis='x',which='minor',direction='in',bottom=True,left=True,right=True,top=True)
    plt.tick_params( axis='y',which='both',direction='in')
    #plt.xticks(np.arange(0.,1.01,0.1))
    plt.ylabel(label_res,fontsize=fos*0.75)
    #PLOT DATA
    for j in range(0,len(tdata)):
      #
      plt.errorbar(tdata[j],res[j],ejdata[j],fmt='.',\
      alpha=0.9,color='#C0C0C0',markersize=rv_markersize,fillstyle='none',zorder=1)
      #
      plt.errorbar(tdata[j],res[j],edata[j],label=telescopes_labels[labvec[j]],\
      fmt=mark[labvec[j]],alpha=1.0,color=rv_colors[labvec[j]], \
      markersize=rv_markersize,fillstyle=rv_fillstyle,zorder=3)
    plt.plot([min(tmodel),max(tmodel)],[0.,0.],'k--',linewidth=1.0,zorder=2)
    #
    yylims = ax1.get_ylim()
    miy = int(max(abs(yylims[0]),abs(yylims[1])))
    if ( miy == 0): #To avoid errors for really small RV variations
      miy = max(abs(yylims[0]),abs(yylims[1]))
    plt.yticks(np.arange(-miy,miy,2.*miy/4.))
    #plt.ylim(-miy,miy)
    plt.minorticks_on()
    plt.xlim(min(tmodel),max(tmodel))
  #
  plt.savefig(fname,format='pdf',bbox_inches='tight')
  plt.savefig(fname[:-3]+'png',format='png',bbox_inches='tight',dpi=300)
  plt.close()

#===========================================================
#        Parameters to be used in plots
#===========================================================

u1_val = np.ndarray(nbands)
u2_val = np.ndarray(nbands)
my_ldc = []
for o in range(0,nbands):
  u1_val[o] =best_value(u1_vec[o],maxloglike,get_value)
  u2_val[o] =best_value(u2_vec[o],maxloglike,get_value)
  my_ldc.append([u1_val[o],u2_val[o]])
my_ldc = np.concatenate(my_ldc)

flag = [False]*4

v_val = [None]*nt
#3 + npars + ldc
for o in range(0,nt):
  v_val[o] = best_value(rv_vec[o],maxloglike,get_value)
  if ( is_log_rv0 ):
    v_val[o] = 10.0**(v_val[o])

t0_val = np.ndarray(nplanets)
P_val  = np.ndarray(nplanets)
e_val  = np.ndarray(nplanets)
w_val  = np.ndarray(nplanets)
i_val  = np.ndarray(nplanets)
a_val  = np.ndarray(nplanets)
tp_val = np.ndarray(nplanets)
k_val  = np.ndarray(nplanets)
rp_val = np.ndarray((nplanets*nradius))

for m in range(0,nplanets):
  t0_val[m] = best_value(T0_vec[m],maxloglike,get_value)
  P_val[m]  = best_value(P_vec[m],maxloglike,get_value)
  e_val[m]  = best_value(e_vec[m],maxloglike,get_value)
  w_val[m]  = best_value(w_vec[m],maxloglike,get_value)
  i_val[m]  = best_value(i_vec[m],maxloglike,get_value)
  a_val[m]  = best_value(ar_vec[m],maxloglike,get_value)
  k_val[m]  = best_value(k_vec[m],maxloglike,get_value)
  tp_val[m] = pti.find_tp(t0_val[m],e_val[m],w_val[m],P_val[m])
  for o in range(0,nradius):
    rp_val[m*nradius+o] = best_value(rr_vec[m*nradius+o],maxloglike,get_value)

alpha_val = best_value(params[4+strends],maxloglike,get_value)
beta_val  = best_value(params[4+strends+1],maxloglike,get_value)

jrv = [None]*n_jrv
for o in range(0,n_jrv):
    jrv[o] = best_value(params[4+sjitrv+o],maxloglike,get_value)

jtr = [None]*n_jtr
for o in range(0,n_jtr):
    jtr[o] = best_value(params[4+sjittr+o],maxloglike,get_value)

pk_rv = []
for m in range(0,np_rv):
    pk_rv.append(best_value(params[4+skrv+m],maxloglike,get_value))

pk_tr = []
for m in range(0,np_tr):
    pk_tr.append(best_value(params[4+sktr+m],maxloglike,get_value))
#===========================================================
#                   Histogram plots
#===========================================================

#Define the labels to be used in the plots

labs = []
elab = '$e$'
wlab = '$\omega$'
ilab = '$i$ (deg)'
klab = '$K$'
alab = '$a/R_\star$'
if ( is_ew ):
    elab = '$\sqrt{e} \sin \omega$'
    wlab = '$\sqrt{e} \cos \omega$'
if ( is_b_factor ):
    ilab = 'b'
if ( is_log_k ):
    klab = '$\log_{10} K$'
if ( is_den_a ):
    alab = '$\\rho_{\star}^{-3}$'
#planet parameter labels
for o in range(0,nplanets):
  etiquetas = ['$T0$'+plabels[o],'$P$'+plabels[o],elab+plabels[o], \
               wlab+plabels[o],ilab+plabels[o],alab+plabels[o], \
               klab+plabels[o]]
  labs.append(etiquetas)
#planet radius labels
for o in range(0,nplanets):
  for m in range(0,nradius):
    labs.append(['$R_p/R_\star$'+plabels[o]+bands[m]])
#LDC labels
for m in range(0,nbands):
  labs.append(['$q_1$'+bands[m],'$q_2$'+bands[m]])
#RV instrument labels
labs.append(telescopes_labels)
#jitter labels
for o in range(0,n_jrv):
  labs.append(['RV_jitter'+str(telescopes_labels[o])])
for o in range(0,n_jtr):
  labs.append(['TR_jitter'+str(bands[o])])
#trends labels
labs.append(['Linear trend'])
labs.append(['Quadratic trend'])
labs.append(krv_labels)
labs.append(ktr_labels)
#Total labels vector
labels = np.concatenate(labs)

def plot_correlations():
  create_plot_correlation(params[4:],labels,col='blue',num=plot_parameters)

def plot_posterior():
    create_plot_posterior(params[4:],labels, cbars='red', nb=50, num=plot_parameters)

def create_plot_posterior(params,plabs,cbars='red',nb=50,num=[]):
  if ( len(num) < 2 ):
    n = range(0,len(params))
  else:
    n = num

  #priorf = np.concatenate([fit_all,fit_ldc,fit_rvs])
  #priorl = np.concatenate([limits,limits_ldc,limits_rvs])
  priorf = prior_flags
  priorl = prior_vals

  plt.figure(1,figsize=(12,4*(len(n))/n_columns_posterior))
  gs = gridspec.GridSpec(nrows=int((len(n)+n_columns_posterior-1)/n_columns_posterior),ncols=n_columns_posterior)
  gs.update(wspace=0.025)
  j = int(0)
  for i in n:
    ax0 = plt.subplot(gs[j])
    vpar, lpar, rpar = find_vals_perc(params[i],1.0)
    moda = my_mode(params[i])
    #best_val = params[i][minchi2_index]
    #plt.axvline(x=best_val,c='yellow')
    plt.axvline(x=vpar,c=cbars)
    plt.axvline(x=moda,c='y',ls='-.')
    plt.axvline(x=vpar-lpar,c=cbars,ls='--')
    plt.axvline(x=vpar+rpar,c=cbars,ls='--')
    plt.xlabel(plabs[i])
    if ( j % n_columns_posterior == 0 ): plt.ylabel('Frequency')
    plt.tick_params( axis='y',which='both',direction='in',labelleft=False)
    plt.tick_params( axis='x',which='both',direction='in')
    if ( is_seaborn_plot ):
      #sns.kdeplot(params[i],label='Posterior, P(M|D)')
      plt.hist(params[i],density=True,bins=nb,histtype='step',label='Posterior, P(M|D)')
    else:
      plt.hist(params[i],density=True,bins=nb,histtype='step',label='Posterior, P(M|D)')
    #Let us plot the prior ranges over the posterior distributions
    if is_plot_prior:
      lx,rx = ax0.get_xlim()
      if (priorf[i] == 'm'  and lx < 0 ): lx = 1e-20 #avoid negative values
      locx = np.arange(lx,rx,(rx-lx)/1000.)
      lp = [None]*len(locx)
      for k in range(0,len(locx)):
        lp[k] = pti.get_priors(priorf[i],[priorl[i*2],priorl[i*2+1]],locx[k])
      plt.plot(locx,lp,alpha=0.8,color='g',label='prior, P(M)')
      #plt.fill_between(locx,lp,alpha=0.3,color='g',label='P(M)')
    #
    if ( i == n[0] ): plt.legend(loc=0, ncol=1,scatterpoints=1,numpoints=1,frameon=True,fontsize=fos*0.5)
    j = int(j + 1)

  fname = outdir+'/'+tic+'_posterior.pdf'
  print('Creating ', fname)
  plt.savefig(fname,format='pdf',bbox_inches='tight')
  fname = outdir+'/'+tic+'_posterior.png'
  plt.savefig(fname,format='png',bbox_inches='tight',dpi=300)
  plt.close()


def create_plot_correlation(params,plabs,col='red',mark='.',num=[]):
  if ( len(num) < 1 ):
    n = range(0,len(params))
  else:
    n = num
  plt.figure(1,figsize=(4*len(n),4*len(n)))
  gs = gridspec.GridSpec(nrows=len(n),ncols=len(n))
  o = 0
  for i in n:
    p = 0
    for j in n:
      if ( j < i ):
        plt.subplot(gs[o*len(n)+p])
        plt.tick_params( axis='y',which='both',direction='in',labelleft=False)
        plt.tick_params( axis='x',which='both',direction='in',labelbottom=False)
        plt.ticklabel_format(useOffset=False, axis='both')
        if ( j == n[0] ):
           plt.ylabel(plabs[i],fontsize=25)
        elif ( j == i - 1 ):
          plt.tick_params( axis='y',which='both',direction='in',labelleft=False)
          plt.tick_params( axis='x',which='both',direction='in',labelbottom=False)
        else:
          plt.tick_params( axis='y',which='both',direction='in',labelleft=False)
          plt.tick_params( axis='x',which='both',direction='in',labelbottom=False)
        if ( i == n[len(n)-1]):
          plt.xlabel(plabs[j],fontsize=25)
        else:
          plt.tick_params( axis='y',which='both',direction='in',labelleft=False)
          plt.tick_params( axis='x',which='both',direction='in',labelbottom=False)
        plt.hist2d(params[j],params[i],bins=100,norm=LogNorm())
        p = p + 1
    o = o + 1

  fname = outdir+'/'+tic+'_correlations.pdf'
  print('Creating ', fname)
  plt.savefig(fname,format='pdf',bbox_inches='tight')
  fname = outdir+'/'+tic+'_correlations.png'
  plt.savefig(fname,format='png',bbox_inches='tight',dpi=300)
  plt.close()

def create_corner_plot():
  import corner

  #update plot_parameters vector
  npp = list(plot_parameters)
  for o in range(0,len(plot_parameters)):
      npp[o] = 4 + plot_parameters[o]

  #Let us take only the values to be plotted
  newpars = [0.0]*len(npp)
  newlabs = [0.0]*len(npp)
  true_params = [0.0]*len(npp)
  for o in range(0,len(npp)):
      newpars[o] = params[npp[o]]
      newlabs[o] = labels[plot_parameters[o]]
      true_params[o] = best_value(newpars[o],maxloglike,get_value)


  #Let us prepare the vector for corner
  data = np.zeros(shape=(len(newpars[0]),len(npp)))
  for o in range(0,len(newpars[0])):
      dumvec = []
      for m in range(0,len(npp)):
          dumvec.append(newpars[m][o])
      data[o] = dumvec

  figure = corner.corner(data, labels=newlabs, \
                       quantiles=[0.16, 0.5, 0.84], \
                        show_titles=True,  )
  fname = outdir+'/'+tic+'_corner.pdf'
  print('Creating ', fname)
  plt.savefig(fname,format='pdf',bbox_inches='tight')
  fname = outdir+'/'+tic+'_corner.png'
  plt.savefig(fname,format='png',bbox_inches='tight',dpi=300)
  plt.close()
