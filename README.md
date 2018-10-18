# Métodos Numéricos 

O projeto aqui presente prevê pontos posteriores de uma função f(x) dado um ponto (ou mais) e sua respectiva função derivada. Tal previsão pode ser feita utilizando um dos metodos seguintes:
- Método de Euler
- Método de Euler Inverso
- Método de Euler Aprimorado
- Método de Runge-Kutta
- Método de Adams-Bashforth
- Método de Adams-Moulton
- Método de Diferenciação Inversa

## Começando

Para ter acesso a esse conteúdo basta ter o GitHub em instalado sua máquina e inserir os seguintes comandos:

-git init
-git clone https://github.com/luizkof/Metodos-Numericos


Pronto já tem todo projeto em sua máquina. Aproveite!

### Pré Requisitos

Tendo o projeto em sua máquina, bastar ter os sequintes requisitos

- Python 2.7.15 (ou versões superiores) já instalado
- As bibliotecas Matplotlib Numpy e Sympy já instaladas


### Instalando
* [Para instalar o MatPlotLib](https://matplotlib.org/users/installing.html) 
* [Para instalar o Sympy](https://docs.sympy.org/latest/install.html) 

## Executando

Para executar é bem simples, basta escrever o a função que deseja calcula pontos posteriores.
É necessário seguir o seguinte padrão: 
- Metodo Euler,Euler Inverso,Euler Aprimorado e Runge-Kutta: Recebe como entrada os valores y(0),t(0), h, quantidade de passos, a função. E
calcule cada passo do método, onde h é o incremento. Exemplos:

```
    euler 0 0 0.1 20 1-t+4*y
    euler_inverso 0 0 0.1 20 1-t+4*y
    euler_aprimorado 0 0 0.1 20 1-t+4*y
    runge_kutta 0 0 0.1 20 1-t+4*y

```

- Mettodo Adam-Bashforth, Adam-Multon, Fórmula Inversa: Recebe como entrada a lista de valores de y,t(0), h, quantidade de passos, a função,
a ordem (de 2 a 8). E calcule cada passo do método.Também pode obter os valores iniciais por outros metodos. Exemplos:

```
    adam_bashforth 0.0 0.1 0.23 0.402 0.6328 0 0.1 20 1-t+4*y 6
    adam_bashforth_by_euler_inverso 0 0 0.1 20 1-t+4*y 6
    adam_multon 0.0 0.1 0.23 0.402 0.6328 0 0.1 20 1-t+4*y 6
    adam_multon_by_runge_kutta 0 0 0.1 20 1-t+4*y 6
    formula_inversa 0.0 0.1 0.23 0.402 0.6328 0 0.1 20 1-t+4*y 6
    formula_inversa_by_euler 0 0 0.1 20 1-t+4*y 6
    formula_inversa_by_runge_kutta 0 0 0.1 20 1-t+4*y 6
```

Os métodos devem ser escrito exatamente dessa forma( mesmo caractéres).
Em seguida basta executar a seguinte linha de comando no terminal dentro da pasta que está o projeto:
```
    python metodos.py
```

## Construido com:

* [Python](https://www.python.org/) 
* [MatPlot](https://matplotlib.org/) 
* [Sympy](https://www.sympy.org/pt/index.html) 



## Versioning
Usamos o GitHub para controle de versão. Veja mais em https://github.com/luizkof/Metodos-Numericos.

## Authors

* **Luiz Albuquerque** - *Initial work* - [PurpleBooth](https://github.com/luizkof)


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details


