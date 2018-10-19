import matplotlib.pyplot as plt
import numpy as np
import os
from sympy import sympify



incrementoHorizontal = 0.001

arq = open('resultado.txt', 'w')


def calculaPontosbyMetodo(YByMetodos,y,x0,ordem,incremento,expr, texto,metodo):

    t0 = float(x0)
    

    #by_Euler
    if(metodo == 1 or metodo == 6 or metodo == 11):
        if(metodo == 1):
            texto.append('Metodo Adan-Bashforth por Euler\n')
            texto.append('y('+ str(x0) + ") = " + str(y)+"\n" )
        elif(metodo == 6):
            texto.append('Metodo Adan-Multon por Euler\n')
            texto.append('y('+ str(x0) + ") = " + str(y) + "\n")
        elif(metodo == 11):
            texto.append('Metodo Formula Inversa de Diferenciacao por Euler\n')
            texto.append('y('+ str(x0) + ") = " + str(y)+"\n" )



        texto.append("h = " + str(incremento) + '\n')
        cont=1
        texto.append(str(0) + " "+ str(y) +'\n')
        while (cont < ordem-1 ):
            
            y = float(y + incremento*expr.subs([("t", t0) , ("y", float(y))]))
            texto.append(str(cont) + " "+ str(y) +'\n')
            t0 = t0+incremento
            YByMetodos.append(y)
            cont = cont+1

            

    elif(metodo == 2 or metodo == 7 or metodo == 12):
        if(metodo == 2):
            texto.append('Metodo Adan-Bashforth por Euler Inverso\n')
            texto.append('y('+ str(x0) + ") = " + str(y)+'\n' )
        elif(metodo == 7):
            texto.append('Metodo Adan-Multon por Euler Inverso\n')
            texto.append('y('+ str(x0) + ") = " + str(y)+'\n' )
        elif(metodo == 12):
            texto.append('Metodo Formula Inversa de Diferenciacao por Euler Inverso\n')
            texto.append('y('+ str(x0) + ") = " + str(y) +'\n')

        
        texto.append("h = " + str(incremento) + '\n')
        cont=1
        texto.append(str(0) + " "+ str(y) +'\n')
        while (cont < ordem-1 ):
            ybyEuler = float(y + incremento*expr.subs([("t", t0) , ("y", float(y))]))
            y = float(y + incremento*expr.subs([("t", t0 + incremento) , ("y", float(ybyEuler))]))
            texto.append(str(cont) + " "+ str(y) +'\n')
            t0 = t0+incremento
            YByMetodos.append(y)
            cont = cont+1

    elif(metodo == 3  or metodo == 8 or metodo ==13):
        if(metodo == 3):
            texto.append('Metodo Adan-Bashforth por Euler Aprimorado\n')
            texto.append('y('+ str(x0) + ") = " + str(y)+"\n" )
        elif(metodo == 8):
            texto.append('Metodo Adan-Multon por Euler Aprimorado\n')
            texto.append('y('+ str(x0) + ") = " + str(y) +"\n")
        elif(metodo == 13):
            texto.append('Metodo Formula Inversa de Diferenciacao por Euler Aprimorado\n')
            texto.append('y('+ str(x0) + ") = " + str(y)+"\n" )

        texto.append("h = " + str(incremento) + '\n')
        cont=1
        texto.append(str(0) + " "+ str(y) +'\n')
        
        while (cont < ordem-1 ):
            ybyEuler = float(y + incremento*expr.subs([("t", t0) , ("y", float(y))]))            
            y = y + ( expr.subs([("t", t0 + incremento) , ("y", float(ybyEuler))]) +  expr.subs([("t", t0) , ("y", float(y))])) * incremento *0.5
            texto.append(str(cont) + " "+ str(y) +'\n')
            t0 = t0+incremento
            YByMetodos.append(y)
            cont = cont+1

    elif(metodo == 4 or metodo == 9 or metodo ==14):
        if(metodo == 4):
            texto.append('Metodo Adan-Bashforth por Runge-Kuuta\n')
            texto.append('y('+ str(x0) + ") = " + str(y) +"\n")
        elif(metodo == 9):
            texto.append('Metodo Adan-Multon por Runge-Kuuta\n')
            texto.append('y('+ str(x0) + ") = " + str(y) +"\n")
        elif(metodo == 14):
            texto.append('Metodo Formula Inversa de Diferenciacao por Runge-Kutta\n')
            texto.append('y('+ str(x0) + ") = " + str(y) +"\n")

        
        texto.append("h = " + str(incremento) + '\n')
        cont=1
        texto.append(str(0) + " "+ str(y) +'\n')
        while (cont < ordem-1 ):

            k1 = float(expr.subs([("t", t0) , ("y", float(y))]))
            k2 = float(expr.subs([("t", t0 + incremento/2) , ("y", float(y + incremento*k1/2))]))
            k3 = float(expr.subs([("t", t0 + incremento/2) , ("y", float(y + incremento*k2/2))]))
            k4 = float(expr.subs([("t", t0 + incremento) , ("y", float(y + incremento*k3))]))
            
            y = y + ((k1+2*k2+2*k3+k4)*incremento)/6
            texto.append(str(cont) + " "+ str(y) +'\n')
            t0 = t0+incremento
            YByMetodos.append(y)
            cont = cont+1

def metodoAdamMulton(y,x0,incremento, qtdPassos, funcaoDerivada, ordem, metodo,id):
    plt.clf()
    expr = sympify(funcaoDerivada)
    plt.title('AdamMulton')
    plt.grid(True)
    texto = []

    t0 = float(x0)
    
    YByMetodos = [0.0]
    YByMetodos.pop(0)
    YByMetodos.append(y)
    if(metodo==5):
        texto.append('Metodo Adan-Multon\n')
        texto.append('y('+ str(x0) + ") = " + str(y[0])+"\n" )
        texto.append("h = " + str(incremento) + '\n')

        cont=0
        while(cont<len(y)):
            a = y[cont]            
            texto.append(str(cont) + " "+ str(a) +'\n')
            cont = cont+1
        YByMetodos = y
    
    
    calculaPontosbyMetodo(YByMetodos,y,x0,ordem,incremento,expr, texto,metodo)
    #add ponto n+1
    YByMetodos.append(0.0)
        
    t0 = float(x0)
    t0 = t0+(ordem-1)*incremento
    
    contador = ordem-1

    
    
    
    while(contador<= qtdPassos):
        #calcula passo seguinte
        if(metodo == 5 or metodo == 6):
            YByMetodos[-1] = YByMetodos[-2] + expr.subs([("t", t0-incremento) , ("y", float(YByMetodos[-2]))])*incremento
        if(metodo ==7 ):
            #apilando metodo de Euler
            yn1 =YByMetodos[-2] + expr.subs([("t", t0-incremento) , ("y", float(YByMetodos[-2]))])*incremento

            #aplicando metodo de Euler inverso  
            YByMetodos[-1] = YByMetodos[-2] + expr.subs([("t", t0) , ("y", float(yn1))])*incremento
        if(metodo ==8 or metodo == 9) :
            #apilando metodo de Euler
            yn1 =YByMetodos[-2] + expr.subs([("t", t0-incremento) , ("y", float(YByMetodos[-2]))])*incremento

            #aplicando metodo de Euler inverso  
            YByMetodos[-1] = YByMetodos[-2] + ( expr.subs([("t", t0) , ("y", float(yn1))]) +  expr.subs([("t", t0-incremento) , ("y", float(YByMetodos[-2]))]))*incremento*0.5    
            
        if(ordem == 2):
            parte1 = ( (1.0/2.0) * float( expr.subs([("t", t0) , ("y", YByMetodos[1])])    ) )
            parte2 =  ( (1.0/2.0) * float(expr.subs([("t", t0-incremento) , ("y", YByMetodos[0])])) )
            novoValor = YByMetodos[ordem-2] + incremento * ( parte1 + parte2 )   
            plt.plot(t0,novoValor,'go--', linewidth=1, markersize=1)
            YByMetodos.pop(0)
            YByMetodos.pop(-1)
            YByMetodos.append(novoValor)
            YByMetodos.append(0.0)
        elif(ordem == 3):
            parte1 = ( (5.0/12.0) * float( expr.subs([("t", t0) , ("y", YByMetodos[2])])    ) )
            parte2 =  ( (2.0/3.0) * float(expr.subs([("t", t0-incremento) , ("y", YByMetodos[1])])) )
            parte3 = ( (-1.0/12.0) * float(expr.subs([("t", t0-2*incremento) , ("y", YByMetodos[0])])) ) 
            novoValor = YByMetodos[ordem-2] + incremento * ( parte1 +parte2 + parte3)   
            plt.plot(t0,novoValor,'go--', linewidth=1, markersize=1)
            YByMetodos.pop(0)
            YByMetodos.pop(-1)
            YByMetodos.append(novoValor)
            YByMetodos.append(0.0)
        elif(ordem == 4):
            parte1 = ( (3.0/8.0) * float( expr.subs([("t", t0) , ("y", YByMetodos[3])])    ) )
            parte2 =  ( (19.0/24.0) * float(expr.subs([("t", t0-incremento) , ("y", YByMetodos[2])])) )
            parte3 = ( (-5.0/24.0) * float(expr.subs([("t", t0-2*incremento) , ("y", YByMetodos[1])])) ) 
            parte4 = ( (1.0/24.0) * float(expr.subs([("t", t0-3*incremento) , ("y", YByMetodos[0])])) ) 
            novoValor = YByMetodos[ordem-2] + incremento * ( parte1 +parte2  + parte3 +parte4)  
            plt.plot(t0,novoValor,'go--', linewidth=1, markersize=1)
            YByMetodos.pop(0)
            YByMetodos.pop(-1)
            YByMetodos.append(novoValor)
            YByMetodos.append(0.0)
        elif(ordem == 5):
            parte1 = ( (251.0/720.0) * float( expr.subs([("t", t0) , ("y", YByMetodos[4])])    ) )
            parte2 =  ( (323.0/360.0) * float(expr.subs([("t", t0-incremento) , ("y", YByMetodos[3])])) )
            parte3 = ( (-11.0/30.0) * float(expr.subs([("t", t0-2*incremento) , ("y", YByMetodos[2])])) ) 
            parte4 = ( (53.0/360.0) * float(expr.subs([("t", t0-3*incremento) , ("y", YByMetodos[1])])) ) 
            parte5 = ( (-19.0/720.0) * float(expr.subs([("t", t0 - 4*incremento) , ("y", YByMetodos[0])])) )
            novoValor = YByMetodos[ordem-2] + incremento * ( parte1 +parte2  + parte3 + parte4 + parte5 ) 
            plt.plot(t0,novoValor,'go--', linewidth=1, markersize=1)
            YByMetodos.pop(0)
            YByMetodos.pop(-1)
            YByMetodos.append(novoValor)
            YByMetodos.append(0.0)
        
        #k=5
        elif(ordem == 6):
            parte1 = ( (95.0/288.0) * float( expr.subs([("t", t0) , ("y", YByMetodos[5])])    ) )
            parte2 =  ( (1427.0/1440.0) * float(expr.subs([("t", t0-incremento) , ("y", YByMetodos[4])])) )
            parte3 = ( (-133.0/240.0) * float(expr.subs([("t", t0 - 2*incremento) , ("y", YByMetodos[3])])) ) 
            parte4 = ( (241.0/720.0) * float(expr.subs([("t", t0 - 3*incremento) , ("y", YByMetodos[2])])) ) 
            parte5 = ( (-173.0/1440.0) * float(expr.subs([("t", t0 - 4*incremento) , ("y", YByMetodos[1])])) )
            parte6 = ( (3.0/160.0) * float(expr.subs([("t", t0 - 5*incremento) , ("y", YByMetodos[0])])) )
            novoValor = YByMetodos[ordem-2] + incremento * ( parte1 +parte2  + parte3 + parte4 + parte5 + parte6)
            plt.plot(t0,novoValor,'go--', linewidth=1, markersize=1)
            YByMetodos.pop(0)
            YByMetodos.pop(-1)
            YByMetodos.append(novoValor)
            YByMetodos.append(0.0)
            
        
        #k=6
        elif(ordem == 7):
            parte1 = ( (19087.0/60480.0) * float( expr.subs([("t", t0) , ("y", YByMetodos[6])])    ) )
            parte2 =  ( (2713.0/2520.0) * float(expr.subs([("t", t0-incremento) , ("y", YByMetodos[5])])) )
            parte3 = ( (-15487.0/20160.0) * float(expr.subs([("t", t0 - 2*incremento) , ("y", YByMetodos[4])])) ) 
            parte4 = ( (586.0/945.0) * float(expr.subs([("t", t0 - 3*incremento) , ("y", YByMetodos[3])])) ) 
            parte5 = ( (-6737.0/20160.0) * float(expr.subs([("t", t0 - 4*incremento) , ("y", YByMetodos[2])])) )
            parte6 = ( (263.0/2520.0) * float(expr.subs([("t", t0 - 5*incremento) , ("y", YByMetodos[1])])) )
            parte7 = ( (-863.0/60480.0) * float(expr.subs([("t", t0 - 6*incremento) , ("y", YByMetodos[0])])) )
            novoValor = YByMetodos[ordem-2] + incremento * ( parte1 +parte2  + parte3 + parte4 + parte5 + parte6 + parte7)
            plt.plot(t0,novoValor,'go--', linewidth=1, markersize=1)
            YByMetodos.pop(0)
            YByMetodos.pop(-1)
            YByMetodos.append(novoValor)
            YByMetodos.append(0.0)
            
        texto.append(str(contador) + " "+ str(YByMetodos[-2]) +'\n')
        contador = contador +1
        t0 = t0+ incremento


        plt.savefig("imgs/" +str(id) + ".png")

    #escreve no arq
    
    arq.writelines(texto)
    arq.write('\n')
    
    #######################


def metodoInverso(y,x0,incremento, qtdPassos, funcaoDerivada, ordem, metodo,id):
    plt.clf()
    expr = sympify(funcaoDerivada)
    plt.title('MetodoInverso')
    plt.grid(True)
    texto = []

    t0 = float(x0)
    
    YByMetodos = [0.0]
    YByMetodos.pop(0)
    YByMetodos.append(y)
    if(metodo==10):
        texto.append('Metodo Formula Inversa de Diferenciacao\n')
        texto.append('y('+ str(x0) + ") = " + str(y[0])+"\n" )
        texto.append("h = " + str(incremento) + '\n')

        cont=0
        while(cont<len(y)):
            a = y[cont]            
            texto.append(str(cont) + " "+ str(a) +'\n')
            cont = cont+1
        YByMetodos = y
    
    
    calculaPontosbyMetodo(YByMetodos,y,x0,ordem,incremento,expr, texto,metodo)
    #add ponto n+1
    YByMetodos.append(0.0)
        
    t0 = float(x0)
    t0 = t0+(ordem-1)*incremento
    
    contador = ordem-1

    
    while(contador<= qtdPassos):
        #calcula passo seguinte
        #Euler
        if(metodo == 10 or metodo == 11):
            YByMetodos[-1] = YByMetodos[-2] + expr.subs([("t", t0-incremento) , ("y", float(YByMetodos[-2]))])*incremento
        #Euler inverso
        if(metodo ==12 ):
            #apilando metodo de Euler
            yn1 =YByMetodos[-2] + expr.subs([("t", t0-incremento) , ("y", float(YByMetodos[-2]))])*incremento

            #aplicando metodo de Euler inverso  
            YByMetodos[-1] = YByMetodos[-2] + expr.subs([("t", t0) , ("y", float(yn1))])*incremento
            
        
        #Euler aprimorado
        if(metodo ==13 or metodo == 14) :
            #apilando metodo de Euler
            yn1 =YByMetodos[-2] + expr.subs([("t", t0-incremento) , ("y", float(YByMetodos[-2]))])*incremento

            #aplicando metodo de Euler Aprimorado
            YByMetodos[-1] = YByMetodos[-2] + ( expr.subs([("t", t0) , ("y", float(yn1))]) +  expr.subs([("t", t0-incremento) , ("y", float(YByMetodos[-2]))]))*incremento*0.5    
                    
        #k=1
        if(ordem == 2):
            parte1 = ( (1.0/1.0) * float( expr.subs([("t", t0) , ("y", YByMetodos[1])])    ) )*incremento
            parte2 = ( (1.0/1.0) * YByMetodos[0])
            novoValor =  ( parte1 + parte2 )   
            plt.plot(t0,novoValor,'go--', linewidth=1, markersize=1)
            YByMetodos.pop(0)
            YByMetodos.pop(-1)
            YByMetodos.append(novoValor)
            YByMetodos.append(0.0)
        
        #k=2
        elif(ordem == 3):
            parte1 = ( (2.0/3.0) * float( expr.subs([("t", t0) , ("y", YByMetodos[2])])    ) )*incremento
            parte2 =  ( (4.0/3.0) * YByMetodos[1])
            parte3 = ( (-1.0/3.0) * YByMetodos[0])
            novoValor = ( parte1 +parte2 + parte3)   
            plt.plot(t0,novoValor,'go--', linewidth=1, markersize=1)
            YByMetodos.pop(0)
            YByMetodos.pop(-1)
            YByMetodos.append(novoValor)
            YByMetodos.append(0.0)

        #k=3
        elif(ordem == 4):
            parte1 = ( (6.0/11.0) * float( expr.subs([("t", t0) , ("y", YByMetodos[3])])    ) )*incremento
            parte2 =  ( (18.0/11.0) * YByMetodos[2])
            parte3 = ( (-9.0/11.0) * YByMetodos[1])
            parte4 = ( (2.0/11.0) * YByMetodos[0])
            novoValor = ( parte1 +parte2  + parte3 +parte4)  
            plt.plot(t0,novoValor,'go--', linewidth=1, markersize=1)
            YByMetodos.pop(0)
            YByMetodos.pop(-1)
            YByMetodos.append(novoValor)
            YByMetodos.append(0.0)

        #k=4
        elif(ordem == 5):
            parte1 = ( (12.0/25.0) * float( expr.subs([("t", t0) , ("y", YByMetodos[4])])    ) )*incremento
            parte2 = ( (48.0/25.0) * YByMetodos[3])
            parte3 = ( (-36.0/25.0) * YByMetodos[2])
            parte4 = ( (16.0/25.0) * YByMetodos[1])
            parte5 = ( (-3.0/25.0) * YByMetodos[0])
            novoValor =  ( parte1 +parte2  + parte3 + parte4 + parte5 ) 
            plt.plot(t0,novoValor,'go--', linewidth=1, markersize=1)
            YByMetodos.pop(0)
            YByMetodos.pop(-1)
            YByMetodos.append(novoValor)
            YByMetodos.append(0.0)
        
        #k=5
        elif(ordem == 6):
            parte1 = ( (60.0/137.0) * float( expr.subs([("t", t0) , ("y", YByMetodos[5])])    ) )*incremento
            parte2 =  ( (300.0/137.0) * YByMetodos[4])
            parte3 = ( (-300.0/137.0) * YByMetodos[3])
            parte4 = ( (200.0/137.0) * YByMetodos[2])
            parte5 = ( (-75.0/137.0) * YByMetodos[1])
            parte6 = ( (12.0/137.0) * YByMetodos[0])
            novoValor =  ( parte1 +parte2  + parte3 + parte4 + parte5 + parte6)
            plt.plot(t0,novoValor,'go--', linewidth=1, markersize=1)
            YByMetodos.pop(0)
            YByMetodos.pop(-1)
            YByMetodos.append(novoValor)
            YByMetodos.append(0.0)
        
        #k=6
        elif(ordem == 7):
            parte1 = ( (60.0/147.0) * float( expr.subs([("t", t0) , ("y", YByMetodos[6])])    ) )*incremento
            parte2 = ( (360.0/147.0) * YByMetodos[5])
            parte3 = ( (-450.0/147.0) * YByMetodos[4])
            parte4 = ( (400.0/147.0) * YByMetodos[3])
            parte5 = ( (-225.0/147.0) * YByMetodos[2])
            parte6 = ( (72.0/147.0) * YByMetodos[1])
            parte7 = ( (-10.0/147.0) * YByMetodos[0])
            novoValor =  ( parte1 +parte2  + parte3 + parte4 + parte5 + parte6 + parte7)
            plt.plot(t0,novoValor,'go--', linewidth=1, markersize=1)
            YByMetodos.pop(0)
            YByMetodos.pop(-1)
            YByMetodos.append(novoValor)
            YByMetodos.append(0.0)
            
        texto.append(str(contador) + " "+ str(YByMetodos[-2]) +'\n')
        contador = contador +1
        t0 = t0+ incremento


        plt.savefig("imgs/" +str(id) + ".png")


    #escreve no arq
    
    arq.writelines(texto)
    arq.write('\n')
    
    #######################




def metodoAdamBashforth(y,x0,incremento, qtdPassos, funcaoDerivada, ordem, metodo,id):
    plt.clf()
    expr = sympify(funcaoDerivada)
    

    plt.title('AdamBashforth')
    plt.grid(True)
    texto = []

    t0 = float(x0)
    
    YByMetodos = [0.0]
    YByMetodos.pop(0)
    YByMetodos.append(y)
    if(metodo==0):
        texto.append('Metodo Adan-Bashforth\n')
        texto.append('y('+ str(x0) + ") = " + str(y[0])+"\n" )
        texto.append("h = " + str(incremento) + '\n')

        cont=0
        while(cont<len(y)):
            a = y[cont]            
            texto.append(str(cont) + " "+ str(a) +'\n')
            cont = cont+1
        YByMetodos = y
    
    
    calculaPontosbyMetodo(YByMetodos,y,x0,ordem,incremento,expr, texto,metodo)
        
    t0 = float(x0)
    t0 = t0+(ordem-2)*incremento
    contador = ordem-1
    while(contador<= qtdPassos):
        

        if(ordem == 2):
            parte1 = ( (1.0/1.0) * float( expr.subs([("t", t0) , ("y", YByMetodos[0])])    ) )
            novoValor = YByMetodos[ordem-2] + incremento * ( parte1 )   
            plt.plot(t0,novoValor,'go--', linewidth=1, markersize=1)
            YByMetodos.pop(0)
            YByMetodos.append(novoValor)
        elif(ordem == 3):
            parte1 = ( (3.0/2.0) * float( expr.subs([("t", t0) , ("y", YByMetodos[1])])    ) )
            parte2 =  ( (-1.0/2.0) * float(expr.subs([("t", t0-incremento) , ("y", YByMetodos[0])])) )
            novoValor = YByMetodos[ordem-2] + incremento * ( parte1 +parte2 )   
            plt.plot(t0,novoValor,'go--', linewidth=1, markersize=1)
            YByMetodos.pop(0)
            YByMetodos.append(novoValor)
        elif(ordem == 4):
            parte1 = ( (23.0/12.0) * float( expr.subs([("t", t0) , ("y", YByMetodos[2])])    ) )
            parte2 =  ( (-4.0/3.0) * float(expr.subs([("t", t0-incremento) , ("y", YByMetodos[1])])) )
            parte3 = ( (5.0/12.0) * float(expr.subs([("t", t0-2*incremento) , ("y", YByMetodos[0])])) ) 
            novoValor = YByMetodos[ordem-2] + incremento * ( parte1 +parte2  + parte3)  
            plt.plot(t0,novoValor,'go--', linewidth=1, markersize=1)
            YByMetodos.pop(0)
            YByMetodos.append(novoValor)
        elif(ordem == 5):
            parte1 = ( (55.0/24.0) * float( expr.subs([("t", t0) , ("y", YByMetodos[3])])    ) )
            parte2 =  ( (-59.0/24.0) * float(expr.subs([("t", t0-incremento) , ("y", YByMetodos[2])])) )
            parte3 = ( (37.0/24.0) * float(expr.subs([("t", t0-2*incremento) , ("y", YByMetodos[1])])) ) 
            parte4 = ( (-3.0/8.0) * float(expr.subs([("t", t0-3*incremento) , ("y", YByMetodos[0])])) ) 
            novoValor = YByMetodos[ordem-2] + incremento * ( parte1 +parte2  + parte3 + parte4) 
            plt.plot(t0,novoValor,'go--', linewidth=1, markersize=1)
            YByMetodos.pop(0)
            YByMetodos.append(novoValor)
        elif(ordem == 6):
            parte1 = ( (1901.0/720.0) * float( expr.subs([("t", t0) , ("y", YByMetodos[4])])    ) )
            parte2 =  ( (-1387.0/360.0) * float(expr.subs([("t", t0-incremento) , ("y", YByMetodos[3])])) )
            parte3 = ( (109.0/30.0) * float(expr.subs([("t", t0 - 2*incremento) , ("y", YByMetodos[2])])) ) 
            parte4 = ( (-637.0/360.0) * float(expr.subs([("t", t0 - 3*incremento) , ("y", YByMetodos[1])])) ) 
            parte5 = ( (251.0/720.0) * float(expr.subs([("t", t0 - 4*incremento) , ("y", YByMetodos[0])])) )
            novoValor = YByMetodos[ordem-2] + incremento * ( parte1 +parte2  + parte3 + parte4 + parte5)
            plt.plot(t0,novoValor,'go--', linewidth=1, markersize=1)
            YByMetodos.pop(0)
            YByMetodos.append(novoValor)
            
        elif(ordem == 7):
            parte1 = ( (4277.0/1440.0) * float( expr.subs([("t", t0) , ("y", YByMetodos[5])])    ) )
            parte2 =  ( (-2641.0/480.0) * float(expr.subs([("t", t0-incremento) , ("y", YByMetodos[4])])) )
            parte3 = ( (4991.0/720.0) * float(expr.subs([("t", t0 - 2*incremento) , ("y", YByMetodos[3])])) ) 
            parte4 = ( (-3649.0/720.0) * float(expr.subs([("t", t0 - 3*incremento) , ("y", YByMetodos[2])])) ) 
            parte5 = ( (959.0/480.0) * float(expr.subs([("t", t0 - 4*incremento) , ("y", YByMetodos[1])])) )
            parte6 = ( (-95.0/288.0) * float(expr.subs([("t", t0 - 5*incremento) , ("y", YByMetodos[0])])) )
            novoValor = YByMetodos[ordem-2] + incremento * ( parte1 +parte2  + parte3 + parte4 + parte5)
            plt.plot(t0,novoValor,'go--', linewidth=1, markersize=1)
            YByMetodos.pop(0)
            YByMetodos.append(novoValor)
            
        texto.append(str(contador) + " "+ str(YByMetodos[ordem-2]) +'\n')
        contador = contador +1
        t0 = t0+ incremento


        plt.savefig("imgs/" +str(id) + ".png")

    #escreve no arq
    
    arq.writelines(texto)
    arq.write('\n')
    
    #######################



def metodoRungeKutta(y0,x0,incremento, qtdPassos, funcaoDerivada,id):
    plt.clf()
    expr = sympify(funcaoDerivada)

    plt.title('Runge-Kuuta')
    plt.grid(True)
    texto = []
    texto.append('Metodo de Runge-Kuuta\n')
    texto.append('y('+ str(x0) + ') = ' + str(y0)+ '\n')
    texto.append("h = " + str(incremento) + '\n')

    t0 = float(x0)

    contador =0
    while(contador<=qtdPassos):
        

        k1 = float(expr.subs([("t", t0) , ("y", float(y0))]))
        k2 = float(expr.subs([("t", t0 + incremento/2) , ("y", float(y0 + incremento*k1/2))]))
        k3 = float(expr.subs([("t", t0 + incremento/2) , ("y", float(y0 + incremento*k2/2))]))
        k4 = float(expr.subs([("t", t0 + incremento) , ("y", float(y0 + incremento*k3))]))
        
        texto.append(str(contador) + " "+ str(y0) +'\n')


        y0 = y0 + ((k1+2*k2+2*k3+k4)*incremento)/6
        plt.plot(t0,y0,'go--', linewidth=1, markersize=1)



        contador = contador +1
        t0 = t0+ incremento


    plt.savefig("imgs/" +str(id) + ".png")

    #escreve no arq
    
    arq.writelines(texto)
    arq.write('\n')
    
    #######################






def metodoEuler(y0,x0,incremento, qtdPassos, funcaoDerivada,id):
    
    expr = sympify(funcaoDerivada)

    plt.clf()
    plt.title('Euler')
    plt.grid(True)
    texto = []
    texto.append('Metodo de Euler\n')
    texto.append('y('+ str(x0) + ') = ' + str(y0)+ '\n')
    texto.append("h = " + str(incremento) + '\n')


    t0 = float(x0)

    contador =0
    while(contador<=qtdPassos):
        texto.append(str(contador) + " "+ str(y0) +'\n')
        y0 = y0 + expr.subs([("t", t0) , ("y", float(y0))])*incremento
        plt.plot(t0,y0,'go--', linewidth=1, markersize=1)
        contador = contador +1
        t0 = t0+ incremento

    plt.savefig("imgs/" +str(id) + ".png")

    #escreve no arq
    
    arq.writelines(texto)
    arq.write('\n')
    
    #######################


def metodoEulerInverso(y0,x0,incremento, qtdPassos, funcaoDerivada,id):
    plt.clf()
    expr = sympify(funcaoDerivada)
    plt.title('EulerInverso')
    plt.grid(True)
    texto = []
    texto.append('Metodo de Euler Inverso\n')
    texto.append('y('+ str(x0) + ') = ' + str(y0)+ '\n')
    texto.append("h = " + str(incremento) + '\n')

    t0 = float(x0)

    contador =0
    while(contador<=qtdPassos):

        texto.append(str(contador) + " "+ str(y0) +'\n')

        #apilando metodo de Euler
        yn1 =y0 + expr.subs([("t", t0) , ("y", float(y0))])*incremento

        #aplicando metodo de Euler inverso  
        y0 = y0 + expr.subs([("t", t0+incremento) , ("y", float(yn1))])*incremento
        

        plt.plot(t0,y0,'go--', linewidth=1, markersize=1)

        contador = contador +1
        t0 = t0 + incremento


    plt.savefig("imgs/" +str(id) + ".png")

    #escreve no arq
    
    arq.writelines(texto)
    arq.write('\n')
    
    #######################


def metodoEulerAprimorado(y0,x0,incremento, qtdPassos, funcaoDerivada,id):
    #limpa grafico
    plt.clf()

    expr = sympify(funcaoDerivada)
    plt.title('EulerAprimorado')
    plt.grid(True)
    texto = []
    texto.append('Metodo de Euler Aprimorado\n')
    texto.append('y('+ str(x0) + ') = ' + str(y0)+ '\n')
    texto.append("h = " + str(incremento) + '\n')


    t0 = float(x0)

    contador =0
    while(contador<=qtdPassos):
        
        texto.append(str(contador) + " "+ str(y0) +'\n')
        #apilando metodo de Euler
        yn1 = y0 + expr.subs([("t", t0) , ("y", float(y0))])*incremento

        #aplicando metodo de Euler inverso  
        y0 = y0 + ( expr.subs([("t", t0 + incremento) , ("y", float(yn1))]) +  expr.subs([("t", t0) , ("y", float(y0))]))*incremento*0.5

        plt.plot(t0,y0,'go--', linewidth=1, markersize=1)

        contador = contador + 1
        t0 = t0 + incremento


    plt.savefig("imgs/" +str(id) + ".png")

    #escreve no arq
    
    arq.writelines(texto)
    arq.write('\n')
    
    #######################


def main():

    
    if not(os.path.isdir("./imgs")): # vemos de este diretorio ja existe
        os.mkdir('./imgs')    
    
    
    numeroLinha = 0

    
    try:
        with open('entrada.txt') as arq:
            for linha in arq:
                numeroLinha = numeroLinha +1
                entrada = linha.split()
                if(entrada[0] == 'euler'):
                    metodoEuler(float(entrada[1]), float(entrada[2]), float(entrada[3]), int(entrada[4]), entrada[5],numeroLinha)
                elif(entrada[0] == 'euler_inverso'):
                    metodoEulerInverso(float(entrada[1]), float(entrada[2]), float(entrada[3]), int(entrada[4]), entrada[5],numeroLinha)
                elif(entrada[0] == 'euler_aprimorado'):
                    metodoEulerAprimorado(float(entrada[1]), float(entrada[2]), float(entrada[3]), int(entrada[4]), entrada[5],numeroLinha)
                elif(entrada[0] == 'runge_kutta'):
                    metodoRungeKutta(float(entrada[1]), float(entrada[2]), float(entrada[3]), int(entrada[4]), entrada[5],numeroLinha)
                
                
                
                elif(entrada[0] == 'adam_bashforth'): 
                    listaY = [0.0]
                    listaY.pop(0)
                    cont=1
                    while cont < int(entrada[-1]):
                        listaY.append(float(entrada[cont]))
                        cont = cont+1
                    metodoAdamBashforth(listaY, float(entrada[-5]), float(entrada[-4]), int(entrada[-3]), entrada[-2],int(entrada[-1]),0,numeroLinha)

                elif(entrada[0] == 'adam_bashforth_by_euler'):
                    metodoAdamBashforth(float(entrada[1]), float(entrada[2]), float(entrada[3]), int(entrada[4]), entrada[5],int(entrada[6]),1,numeroLinha)
                elif(entrada[0] == 'adam_bashforth_by_euler_inverso'):
                    metodoAdamBashforth(float(entrada[1]), float(entrada[2]), float(entrada[3]), int(entrada[4]), entrada[5],int(entrada[6]),2,numeroLinha)
                elif(entrada[0] == 'adam_bashforth_by_euler_aprimorado'):
                    metodoAdamBashforth(float(entrada[1]), float(entrada[2]), float(entrada[3]), int(entrada[4]), entrada[5],int(entrada[6]),3,numeroLinha)
                elif(entrada[0] == 'adam_bashforth_by_runge_kutta'):
                    metodoAdamBashforth(float(entrada[1]), float(entrada[2]), float(entrada[3]), int(entrada[4]), entrada[5],int(entrada[6]),4,numeroLinha)
                
                elif(entrada[0] == 'adam_multon'): 
                    listaY = [0.0]
                    listaY.pop(0)
                    cont=1
                    while cont < int(entrada[-1]):
                        listaY.append(float(entrada[cont]))
                        cont = cont+1
                    metodoAdamMulton(listaY, float(entrada[-5]), float(entrada[-4]), int(entrada[-3]), entrada[-2],int(entrada[-1]),5,numeroLinha)

                elif(entrada[0] == 'adam_multon_by_euler'):
                    metodoAdamMulton(float(entrada[1]), float(entrada[2]), float(entrada[3]), int(entrada[4]), entrada[5],int(entrada[6]),6,numeroLinha)
                elif(entrada[0] == 'adam_multon_by_euler_inverso'):
                    metodoAdamMulton(float(entrada[1]), float(entrada[2]), float(entrada[3]), int(entrada[4]), entrada[5],int(entrada[6]),7,numeroLinha)
                elif(entrada[0] == 'adam_multon_by_euler_aprimorado'):
                    metodoAdamMulton(float(entrada[1]), float(entrada[2]), float(entrada[3]), int(entrada[4]), entrada[5],int(entrada[6]),8,numeroLinha)
                elif(entrada[0] == 'adam_multon_by_runge_kutta'):
                    metodoAdamMulton(float(entrada[1]), float(entrada[2]), float(entrada[3]), int(entrada[4]), entrada[5],int(entrada[6]),9,numeroLinha)
                
                elif(entrada[0] == 'formula_inversa'): 
                    listaY = [0.0]
                    listaY.pop(0)
                    cont=1
                    while cont < int(entrada[-1]):
                        listaY.append(float(entrada[cont]))
                        cont = cont+1
                    metodoInverso(listaY, float(entrada[-5]), float(entrada[-4]), int(entrada[-3]), entrada[-2],int(entrada[-1]),10,numeroLinha)

                elif(entrada[0] == 'formula_inversa_by_euler'):
                    metodoInverso(float(entrada[1]), float(entrada[2]), float(entrada[3]), int(entrada[4]), entrada[5],int(entrada[6]),11,numeroLinha)
                elif(entrada[0] == 'formula_inversa_by_euler_inverso'):
                    metodoInverso(float(entrada[1]), float(entrada[2]), float(entrada[3]), int(entrada[4]), entrada[5],int(entrada[6]),12,numeroLinha)
                elif(entrada[0] == 'formula_inversa_by_euler_aprimorado'):
                    metodoInverso(float(entrada[1]), float(entrada[2]), float(entrada[3]), int(entrada[4]), entrada[5],int(entrada[6]),13,numeroLinha)
                elif(entrada[0] == 'formula_inversa_by_runge_kutta'):
                    metodoInverso(float(entrada[1]), float(entrada[2]), float(entrada[3]), int(entrada[4]), entrada[5],int(entrada[6]),14,numeroLinha)
                else:
                    print("Nao existe o metodo " + str(entrada[0]) ) 
                print("Concluiu linha " + str(linha))
        arq.close()
    except IOError:
       print 'Arquivo entrada.txt nao existe'

    print("Confira seu resultado no arquivo resultado.txt")
    
main()
