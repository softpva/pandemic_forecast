# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
"""
Created on Tue May  5 17:40:00 2020

@author: 
pedro eduardo souza abreu pedu.abreu@gmail.com
pedro vitor abreu  eng.pva@gmail.com
"""
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt

def run(
        pop: 'população',
        i_ini: 'infectados iniciais',
        r_ini: 'recuperados iniciais',
        r_zero: 'razão de reprodução',
        gama_inv: 'período de recuperação',
        t_max: 'tempo total analisado',
        data_ini: 'data inicial') -> pd.DataFrame:    
    s_ini = pop - i_ini - r_ini
    gama = 1/gama_inv if gama_inv != 0 else 1.
    beta = gama * r_zero 
    t_list = np.arange(t_max)
    seq_dias = np.arange(data_ini, data_ini + t_max, dtype= 'datetime64[D]')
    res_integr = integrar(s_ini, i_ini, r_ini, beta, gama, t_max)
    s,i,r,i_acum = res_integr
    return pd.DataFrame({'S': s, 'I': i, 'R': r, 'I_acum':i_acum,
                         'D': seq_dias}, index=t_list)

def integrar(
    si: 'suscetíveis iniciais',
    ii: 'infectados iniciais',
    ri: 'resolvendo inicias',
    b: 'taxa de propagação S->I',
    g: 'taxa de recuperação I->R',
    t_max: 'número de dias analisados'
    ): # retorna (s,i,r,i_acum) np.array
    s = np.zeros(t_max, dtype = float)
    i = s.copy(); r = s.copy(); i_acum = s.copy()   
    s[0] = si
    i_acum[0] = ii
    i[0] = ii
    r[0] = 1
    n = si + ii + ri
    for t in np.arange(t_max - 1, dtype = int):
        s[t+1]      = s[t] - (b * s[t] * i[t] / n)
        i[t+1]      = i[t] + (b * s[t] * i[t] / n) - (g * i[t])
        r[t+1]      = r[t] + g * i[t]
        i_acum[t+1] = i_acum[t] + (b * s[t] * i[t] / n)
    return (s,i,r,i_acum)   

def plotar(
        tit: 'Título da figura de plotagem, incui parâmetros principai',
        res: 'DataFrame com resultados',
        data_plt_ini: 'ano-mes-dia inicio da plotagem',
        data_plt_fim: 'ano-mes-dia fim da plotagem'):
    
    fig = plt.figure(edgecolor='b', facecolor='whitesmoke', figsize=(15,7.5))
    ax1 = fig.add_subplot(111)    
    showdate_ini = np.datetime64(data_plt_ini)
    showdate_end = np.datetime64(data_plt_fim)    
    ax1.plot(res['D'], res['I_acum'], 'r',lw=2, label="Infectados Acumuladas")    
    ax1.plot(res['D'], res['I'] , 'r--', lw=2, label="Infectados")
    ax1.plot(res['D'], res['S'], 'b', lw=2, label="Suscetíveis")
    ax1.plot(res['D'], res['R'], 'g', lw=2, label="Recuperados")    
    ax1.set_ylabel('População')
    ax1.set_xlim(showdate_ini,showdate_end)    
    ax1.set_xticks(np.arange(showdate_ini,showdate_end , 7, dtype='datetime64[D]'))
    ax1.set_yticks(np.arange(0, pop, pop/20))
    ax1.xaxis.set_tick_params(rotation=90) 
    ax1.grid(b = True, which = 'major', c = 'silver', lw=1, ls='-') 
    ax1.set_title(tit, fontsize=12, pad =25)        
    legend_1 = ax1.legend(loc= 'center left')
    legend_1.get_frame().set_alpha(0.5)

if __name__ == '__main__':
    pop = 1_000_000
    i_ini,r_ini = 1., 0
    r_zero, gama_inv = 3.,5   
    t_max = 180
    data_ini = np.datetime64('2020-03-01')    
    res = run(pop, i_ini, r_ini, r_zero, gama_inv, t_max, data_ini)    
    tit = 'Fig_1: SIR Padrão Cluster Único. População: {:,} Infect.Ini: {}\
    R0: {} Gama_inv: {}'.format(pop, i_ini, r_zero, gama_inv)    
    plotar(tit,res,'2020-03-01','2020-08-30')