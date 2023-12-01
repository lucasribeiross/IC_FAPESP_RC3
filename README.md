# IC_FAPESP_RC3_codes
Repositório dos códigos utilizados para a apresentação dos resultados do RC3 - IC FAPESP

## Evolução temporal do sistema
Inicialmente estamos interessados em saber qual o tempo de dinâmica do sistema é necessário para que o estado estacionário seja alcançado. Podemos então calcular a dinâmica do sistema e verificar a região no tempo para obter o estado estacionário independente do $\Omega_C$ adotado.

Código: time_evolution.py

## Espectro de transmissão do EIT
Vizualizamos a transmissão do sistema pela dessintonia $\Delta_P$ fixando a frequência $\Omega_C=\kappa$ e analisando como o acoplamento ($g_0$) influencia na posição dos estados vestidos. Em especial, os próximos resultados utilizarão a construção desse espectro.

Código: EIT_transmission.py

## FWHM do espectro de transmissão
Calculando para um regime de frequências de $\Omega_C$ o espectro de transmissão do EIT. São obtidos poucos pontos em torno da largura de meia-altura (FWHM), utilizando uma interpolação polinomial (Spline-SciPy), a fim de otimizar o cálculo numérico. O código se encontra generalizado para $N_{\text{at}}$ número de átomos. A base de Fock $N=6$ foi escolhida truncando as probabilidades do estado de fótons que entram na cavidade segundo a força do campo de bombeio $\epsilon_0$.

Código: test.py
