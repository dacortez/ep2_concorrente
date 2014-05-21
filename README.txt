========================================
EP2 - MAC0438 Programaçao Concorrente

Daniel Augusto Cortez

22/05/2014
========================================


COMPILAÇÃO
----------

Basta utilizar o Makefile

	$> make ep2

O programa utiliza a bibliotema GMP para aritmética de múltipla precisão.
A mesma deve estar copilada e disponível no local padrão.

O programa também utiliza POSIX pthreads e uma barreira de sincronização lá
disponível.


UTILIZAÇÃO
----------

	$> ./ep2 <threads> <f|m> <eps> [d|s|x|y] [p]

Onde:

  threads 	 número de threads (0 utiliza o número de núcles)
        f 	 para com diferença menor do que eps
        m 	 para com termo menor do que eps (ex: 1e-100)
      eps 	 valor que define a precisão
        d 	 informações de debug
        s 	 execução sequencial
        x 	 experimentos sequencial
        y 	 experimentos paralelo
        p 	 parâmetro de compressão da série de Taylor

O valor da presição eps deve ser informada em notação cientĩfica, por exemplo,
1e-100. 

Para detalhes sobre a implementação, consulte o relatório.pdf.
