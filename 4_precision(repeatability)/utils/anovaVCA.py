import pandas as pd
from statsmodels.graphics.factorplots import interaction_plot
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

class VCA:
    def __init__(self,data):
        self.datafile = data     
        if 'xlsx' in self.datafile:
            self.dataframe = pd.read_excel(self.datafile)
        elif 'csv' in self.datafile:
            self.dataframe = pd.read_csv(self.datafile)
        
    def nested_ANOVA(self,arg="A~Day/Run",save_xlsx=False,repeatability_precision=False,brief=True):
        arg_y = arg.split('~')[0]
        arg_gp = arg.split("~")[1].split('/')[0]
        arg_subgp = arg.split("~")[1].split('/')[1]
        
        if 'xlsx' in self.datafile:
            data = pd.read_excel(self.datafile)
        elif 'csv' in self.datafile:
            data = pd.read_csv(self.datafile)
        
        
        group_number = len(data[arg_gp].unique())
        subgroup_number = len(data[arg_subgp].unique())
        a = group_number #天數
        b = subgroup_number  #實驗次數
        n = len(data) / (len(data[arg_gp].unique())*len(data[arg_subgp].unique())) #重複次數
        N = a*b*n


        sample = arg_y

        #SS_group 
        Ti2_bn = 0

        for i in data[arg_gp].unique():
            Ti2_bn += sum(data[data[arg_gp]==i][sample])**2/(b*n)
        G2 = sum(data[sample])**2
        G2_N = G2/N
        SS_group = Ti2_bn - G2_N 


        #SS_subgroup
        Sij2_n = 0
        for i in data[arg_gp].unique():
            for j in data[arg_subgp].unique():
                Sij2_n += sum(data[(data[arg_gp]==i)&(data[arg_subgp]==j)][sample])**2/n
        #print(Sij2_n)
        #print(Ti2_bn)
        SS_subgroup = Sij2_n - Ti2_bn

        #SS_total
        SS_total = sum(data[sample]**2) -  G2_N
        SS_error = SS_total - SS_group - SS_subgroup
        
        #df
        df_gp = a-1
        df_subgp = (b-1)*a
        df_error = (n-1)*a*b
        df_total = N-1
        
        #MS
        MS_gp = SS_group / df_gp
        MS_subgp = SS_subgroup / df_subgp
        MS_error = SS_error / df_error 

        #F
        F_gp = MS_gp / MS_subgp
        F_subgp = MS_subgp / MS_error

        #p-value
        p_gp = stats.f.sf(F_gp, df_gp, df_subgp)
        p_subgp = stats.f.sf(F_subgp, df_subgp, df_error)
        
        #VC
        VC_error = MS_error
        VC_subgp = (MS_subgp-MS_error)/n
        VC_gp = (MS_gp-MS_subgp)/(n*b)
        VC_total = VC_error + VC_subgp + VC_gp
        
        #CV%
        CV_gp = np.sqrt(VC_gp)*100/np.mean(data[sample])
        CV_subgp = np.sqrt(VC_subgp)*100/np.mean(data[sample])
        CV_error = np.sqrt(VC_error)*100/np.mean(data[sample])
        CV_total = np.sqrt(VC_total)*100/np.mean(data[sample])
        
        
        #95%CI of SD
        df_R = N - a*b
        df_WL = (1*MS_gp+1*MS_subgp+2*MS_error)**2/((1*MS_gp)**2/df_gp + (1*MS_subgp)**2/df_subgp + (2*MS_error)**2/df_error)
        SR_lower = np.sqrt(VC_error)*np.sqrt(df_R/stats.chi2.isf(0.05,df_R))
        SR_upper = np.sqrt(VC_error)*np.sqrt(df_R/stats.chi2.isf(0.95,df_R))     
        SWL_lower = np.sqrt(VC_total)*np.sqrt(df_WL/stats.chi2.isf(0.05,df_WL))
        SWL_upper = np.sqrt(VC_total)*np.sqrt(df_WL/stats.chi2.isf(0.95,df_WL))
        Sgp_lower = np.sqrt(VC_gp)*np.sqrt(df_gp/stats.chi2.isf(0.05,df_gp))
        Sgp_upper = np.sqrt(VC_gp)*np.sqrt(df_gp/stats.chi2.isf(0.95,df_gp))
        Ssubgp_lower = np.sqrt(VC_subgp)*np.sqrt(df_subgp/stats.chi2.isf(0.05,df_subgp))
        Ssubgp_upper = np.sqrt(VC_subgp)*np.sqrt(df_subgp/stats.chi2.isf(0.95,df_subgp))

        #95%CI of CV%
        df_R = N - a*b
        df_WL = (1*MS_gp+1*MS_subgp+2*MS_error)**2/((1*MS_gp)**2/df_gp + (1*MS_subgp)**2/df_subgp + (2*MS_error)**2/df_error)
        CVR_lower = SR_lower*100/np.mean(data[sample])
        CVR_upper = SR_upper*100/np.mean(data[sample])
        CVWL_lower = SWL_lower*100/np.mean(data[sample])
        CVWL_upper = SWL_upper*100/np.mean(data[sample])
        CVgp_lower = Sgp_lower*100/np.mean(data[sample])
        CVgp_upper = Sgp_upper*100/np.mean(data[sample])
        CVsubgp_lower = Ssubgp_lower*100/np.mean(data[sample])
        CVsubgp_upper = Ssubgp_upper*100/np.mean(data[sample])
        #用dictionary做Dataframe
        results = {'SS':[SS_group, SS_subgroup, SS_error, SS_total],
                   'MS':[MS_gp,MS_subgp,MS_error,SS_total/df_total],
                   'df':[df_gp, df_subgp, df_error, df_WL],
                   'F':[F_gp, F_subgp,'', ''],
                    'p-value':[p_gp, p_subgp, '', ''],
                   'Variance':[VC_gp,VC_subgp,VC_error,VC_total],
                   'SD':[np.sqrt(VC_gp),np.sqrt(VC_subgp),np.sqrt(VC_error),np.sqrt(VC_total)],
                   'CV%':[CV_gp,CV_subgp,CV_error,CV_total],
                   '95%CI of SD':['%.4f - %.4f'%(Sgp_lower,Sgp_upper),'%.4f - %.4f'%(Ssubgp_lower,Ssubgp_upper),'%.4f - %.4f'%(SR_lower,SR_upper),'%.4f - %.4f'%(SWL_lower,SWL_upper)],
                   '95%CI of CV%':['%.4f - %.4f'%(CVgp_lower,CVgp_upper),'%.4f - %.4f'%(CVsubgp_lower,CVsubgp_upper),'%.4f - %.4f'%(CVR_lower,CVR_upper),'%.4f - %.4f'%(CVWL_lower,CVWL_upper)]}
        columns=['SS','MS', 'df', 'F', 'p-value','Variance','SD','CV%','95%CI of SD','95%CI of CV%']

        results_brief = {'SS':['%.1f'%SS_group, '%.1f'%SS_subgroup, '%.1f'%SS_error, '%.1f'%SS_total],
                   'MS':['%.1f'%MS_gp,'%.1f'%MS_subgp,'%.1f'%MS_error,'%.1f'%(SS_total/df_total)],
                   'df':['%d'%df_gp, '%d'%df_subgp, '%d'%df_error, '%.2f'%df_WL],
                   'F':['%.3f'%F_gp, '%.3f'%F_subgp,'', ''],
                    'p-value':['%.4f'%p_gp, '%.4f'%p_subgp, '', ''],
                   'Variance':['%.3f'%VC_gp,'%.3f'%VC_subgp,'%.3f'%VC_error,'%.3f'%VC_total],
                   'SD':['%.3f'%np.sqrt(VC_gp),'%.3f'%np.sqrt(VC_subgp),'%.3f'%np.sqrt(VC_error),'%.3f'%np.sqrt(VC_total)],
                   'CV%':['%.3f'%CV_gp,'%.3f'%CV_subgp,'%.3f'%CV_error,'%.3f'%CV_total],
                   '95%CI of SD':['%.4f - %.4f'%(Sgp_lower,Sgp_upper),'%.4f - %.4f'%(Ssubgp_lower,Ssubgp_upper),'%.4f - %.4f'%(SR_lower,SR_upper),'%.4f - %.4f'%(SWL_lower,SWL_upper)],
                   '95%CI of CV%':['%.4f - %.4f'%(CVgp_lower,CVgp_upper),'%.4f - %.4f'%(CVsubgp_lower,CVsubgp_upper),'%.4f - %.4f'%(CVR_lower,CVR_upper),'%.4f - %.4f'%(CVWL_lower,CVWL_upper)]}

        aov_table1_brief = pd.DataFrame(results_brief, columns=columns,
                                  index=['%s'%arg_gp, '%s:%s'%(arg_gp,arg_subgp), 'Error','Total'])
        aov_table1 = pd.DataFrame(results, columns=columns,
                                  index=['%s'%arg_gp, '%s:%s'%(arg_gp,arg_subgp), 'Error(repeatability)','Tota(within-laboratory precision)l'])                          
        
        
        
        
        #save ANOVA result?
        if save_xlsx == True:
            aov_table1.to_excel("%s_%s_%s.xlsx"%(arg_y,arg_gp,arg_subgp))
            
        #repeatability_precision
        if repeatability_precision == True:
            S_R = np.sqrt(VC_error)
            S_WL = np.sqrt(VC_total)
            df_R = N - a*b

            print("CV_R%%:%f%%"%(S_R/(sum(data[sample])/N)*100))
            print("CV_WL%%:%f%%"%(S_WL/(sum(data[sample])/N)*100))
            
            
            ########待考證########
            a1 = 1
            a2 = 1
            a3 = 2
            df_WL = (a1*MS_gp+a2*MS_subgp+a3*MS_error)**2/((a1*MS_gp)**2/df_gp + (a2*MS_subgp)**2/df_subgp + (a3*MS_error)**2/df_error)
            ########待考證#########
            
            
            SR_lower = S_R*np.sqrt(df_R/stats.chi2.isf(0.05,df_R))
            SR_upper = S_R*np.sqrt(df_R/stats.chi2.isf(0.95,df_R))
            print("95%%CI of SR:[%.4f-%.4f]"%(SR_lower,SR_upper))
            
            SWL_lower = S_WL*np.sqrt(df_WL/stats.chi2.isf(0.05,df_WL))
            SWL_upper = S_WL*np.sqrt(df_WL/stats.chi2.isf(0.95,df_WL))
            print("95%%CI of SWL:[%.4f-%.4f]"%(SWL_lower,SWL_upper))
        
        if brief == True:
            return aov_table1_brief
        elif brief == False:
            return aov_table1
    
    def plot_meanRun(self,x='Day',y='A',color=True,allplot=False,allcon=['A','B','C','D']):
        if allplot == False:
            if 'xlsx' in self.datafile:
                data = pd.read_excel(self.datafile)
            elif 'csv' in self.datafile:
                data = pd.read_csv(self.datafile)
            days =list(range(1,max(data[x])+1))

            Run1 = []
            Run2 = []
            for i in days:
                Run1.append(np.mean(data[(data[x]==i)&(data['Run']==1)][y]))
                Run2.append(np.mean(data[(data[x]==i)&(data['Run']==2)][y]))


            fig = plt.figure()
            ax1 = fig.add_subplot(1,1,1)
            if color == True:
                plt.plot(Run1)
                plt.plot(Run2) 
            else:
                plt.plot(Run1,color='black')
                plt.plot(Run2,'--',color='black')
            plt.xticks(range(len(days)),days,rotation=0,fontsize='small')
            plt.xlabel('day')
            plt.ylabel('concentration')
            plt.title('substrate %s'%y)
            ax1.legend(['mean of run1','mean of run2'])
            plt.show()  
        
        elif allplot == True:
            fig = plt.figure()
            ax1 = fig.add_subplot(1,1,1) 

            for i in allcon:
                if 'xlsx' in self.datafile:
                    data = pd.read_excel(self.datafile)
                elif 'csv' in self.datafile:
                    data = pd.read_csv(self.datafile)
                days = list(range(1,max(data[x])+1))
                Run1 = []
                Run2 = []    
                for j in days:
                    Run1.append(np.mean(data[(data[x]==j)&(data['Run']==1)][i]))
                    Run2.append(np.mean(data[(data[x]==j)&(data['Run']==2)][i]))

                plt.plot(days,Run1)
                plt.plot(days,Run2)   
            plt.xticks(days,days,rotation=0,fontsize='small')
            plt.xlabel('day')
            plt.ylabel('concentration')
            legendlist = []
            for i in allcon:
                legendlist.append('mean of %s in run1'%i)
                legendlist.append('mean of %s in run2'%i)
            ax1.legend(legendlist,loc='center left', bbox_to_anchor=(1, 0.5))
            plt.show() 
        
    def scatter(self,x='Day',y='A',color=True,allplot=False,allcon=['A','B','C','D']):
        if 'xlsx' in self.datafile:
            data = pd.read_excel(self.datafile)
        elif 'csv' in self.datafile:
            data = pd.read_csv(self.datafile)
        days = list(range(1,max(data[x])+1))

        Run1 = []
        Run2 = []    
        if allplot == False:
            for i in days:
                Run1 += list(data[(data[x]==i)&(data['Run']==1)][y])
                Run2 += list(data[(data[x]==i)&(data['Run']==2)][y])

            days_x = days + days
            days_x.sort()

            fig = plt.figure()
            ax1 = fig.add_subplot(1,1,1)   
            plt.scatter(days_x,Run1)
            plt.scatter(days_x,Run2)
            plt.xticks(days,days,rotation=0,fontsize='small')
            plt.xlabel('day')
            plt.ylabel('concentration')
            plt.title('substrate %s'%y)
            ax1.legend(['Run1','Run2'])
            plt.show() 
            
        elif allplot == True:
            fig = plt.figure()
            ax1 = fig.add_subplot(1,1,1) 
            for i in allcon:
                if 'xlsx' in self.datafile:
                    data = pd.read_excel(self.datafile)
                elif 'csv' in self.datafile:
                    data = pd.read_csv(self.datafile)
                days = list(range(1,max(data[x])+1))
                Run1 = []
                Run2 = []    
                for j in days:
                    Run1 += list(data[(data[x]==j)&(data['Run']==1)][i])
                    Run2 += list(data[(data[x]==j)&(data['Run']==2)][i])
                days_x = days + days
                days_x.sort()
                plt.scatter(days_x,Run1)
                plt.scatter(days_x,Run2)   
            plt.xticks(days,days,rotation=0,fontsize='small')
            plt.xlabel('day')
            plt.ylabel('concentration')
            legendlist = []
            for i in allcon:
                legendlist += ['Run1 of %s'%i,'Run2 of %s'%i]
            ax1.legend(legendlist,loc='center left', bbox_to_anchor=(1, 0.5))
            plt.show() 
