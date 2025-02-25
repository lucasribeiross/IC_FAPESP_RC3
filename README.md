# DD_FAPESP_codes&data
Repositório dos códigos e dados obtidos durante o Doutorado Direto FAPESP Processo 2024/02604-5 - Processamento de Informação Quântica em Eletrodinâmica Quântica de Cavidades.

## Fundamental limit to cavity linewidth narrowing with single atoms

Disponível no arXiv pelo seguinte endereço eletrônico: https://doi.org/10.48550/arXiv.2411.12422.

### Evolução temporal do sistema
Inicialmente estamos interessados em saber qual o tempo de dinâmica do sistema é necessário para que o estado estacionário seja alcançado. Podemos então calcular a dinâmica do sistema e verificar a região no tempo para obter o estado estacionário independente do $\Omega_C$ adotado.

Código modelo quântico: time_evolution.py

Código modelo semi-clássico: semiclassical_EIT_timeevolution.mlx

### Espectro de transmissão do EIT
Vizualizamos a transmissão do sistema pela dessintonia $\Delta_P$ fixando a frequência $\Omega_C=\kappa$ e analisando como o acoplamento ($g_0$) influencia na posição dos estados vestidos. Em especial, os próximos resultados utilizarão a construção desse espectro.

Código modelo quântico: EIT_transmission.py

Código modelo semi-clássico: semiclassical_EIT_spectrumtransmission.mlx

### FWHM do espectro de transmissão
Calculando para um regime de frequências de $\Omega_C$ o espectro de transmissão do EIT. São obtidos poucos pontos em torno da largura de meia-altura (FWHM), utilizando uma interpolação polinomial (Spline-SciPy), a fim de otimizar o cálculo numérico. O código se encontra generalizado para $N_{\text{at}}$ átomos. A base de Fock $N=6$ foi escolhida truncando as probabilidades do estado de fótons que entram na cavidade segundo a força do campo de bombeio $\varepsilon_0$.

Código modelo quântico: test.py

Código modelo semi-clássico: semiclassical_EIT_FWHM.mlx

### Estatística do campo
Alguns dos códigos desenvolvidos para o modelo quântico a fim de analisar a estatística do campo após a interação. O objetivo é identificar uma mudança na natureza do campo decorrente do crescimento de $N_{\text{at}}$. Foi calculada a função de correlação de segunda ordem e as projeções $\langle P_n \rangle$ do projetor correspondente ao número de fótons.

Códigos: Nat_2correlationfunction.py e Pn_Dp.py


## Estimativa do número médio de fótons pelo Algoritmo de Estimação de Fase

