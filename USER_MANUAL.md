# Manual do Utilizador

Este manual explica como usar a aplicação e como interpretar os principais conceitos epidemiológicos e estatísticos apresentados nos separadores. A aplicação é uma ferramenta exploratória não oficial: ajuda a analisar dados de mortalidade do INE, mas não substitui validação epidemiológica, análise clínica, revisão metodológica ou produtos estatísticos oficiais.

Para detalhes técnicos das fórmulas, indicadores e implementação, ver [METHODOLOGY.md](METHODOLOGY.md).

## 1. Visão Geral

A aplicação permite:

- consultar mortalidade observada por ano, local, causa de morte, sexo e grupo populacional;
- comparar Portugal, Norte e uma localização adicional num ano específico;
- calcular óbitos, mortalidade bruta, mortalidade padronizada, mortalidade proporcional e AVPP;
- gerar previsões exploratórias até 30 anos no futuro;
- comparar modelos de previsão;
- avaliar diagnósticos, erros de previsão e possíveis quebras estruturais;
- verificar a disponibilidade de dados nos ficheiros RDS antes de carregar uma análise;
- exportar tabelas em CSV e gráficos em PNG.

A ordem dos separadores é:

1. `Mortalidade Observada`
2. `Previsão Guiada`
3. `Previsão Avançada`
4. `Métricas Anuais`
5. `Disponibilidade de Dados`

## 2. Como Começar

Depois de abrir a aplicação, escolha primeiro a fonte de dados.

`Ficheiros RDS` é a opção recomendada para uso normal. Lê os ficheiros já preparados no repositório ou, se necessário, na localização remota configurada no GitHub. É muito mais rápido do que consultar o INE em directo.

`INE em directo` consulta os indicadores do INE através do pacote `ineptr2`. Deve ser usado quando os ficheiros RDS não contêm os anos, locais ou causas pretendidos, ou quando se quer verificar dados mais recentes. Pode ser bastante lento, sobretudo para o indicador histórico de óbitos `0008206`.

Em geral:

- use `Ficheiros RDS` para análise, exploração e apresentações;
- use `INE em directo` para actualização ou validação pontual;
- consulte `Disponibilidade de Dados` quando não tiver a certeza se os RDS têm a combinação pretendida.

## 3. Conceitos Comuns aos Separadores

### Local de Residência

Define a área geográfica usada na análise. Pode seleccionar uma única localização, como `Alentejo`, ou várias localizações.

Quando selecciona mais do que uma localização, a aplicação soma os óbitos e a população dessas localizações antes de calcular as taxas. Isto é útil para criar uma área agregada, por exemplo uma área de influência ou um conjunto de concelhos.

O campo `Nome da Selecção (opcional)` permite dar um nome a essa agregação. Se for deixado em branco, a aplicação usa uma designação automática.

### Causa de Morte

Define a causa analisada. Algumas análises permitem apenas uma causa por carregamento; outras permitem várias causas.

`Todas as causas de morte` é especialmente importante porque serve como denominador para a `Mortalidade Proporcional`.

### Sexo

Permite analisar:

- homens e mulheres em conjunto;
- apenas homens;
- apenas mulheres.

Quando escolhe ambos os sexos, a aplicação soma os óbitos e a população antes de calcular os indicadores.

### População

`Total` usa todos os grupos etários disponíveis.

`Menos de 75 anos` exclui os grupos etários a partir dos 75 anos. Esta opção pode ser útil quando o foco é mortalidade prematura ou potencialmente evitável, mas deixa de representar a mortalidade total da população.

### Taxa

Nas análises temporais e nas previsões, a taxa seleccionada define a série a observar ou prever. As opções principais são:

- mortalidade bruta;
- mortalidade padronizada;
- mortalidade proporcional, quando aplicável.

## 4. Métricas de Mortalidade

### Óbitos

Número absoluto de mortes para a combinação seleccionada de ano, local, causa, sexo e idade.

Vantagens:

- é a contagem mais directa;
- é útil para planeamento operacional e dimensão do problema.

Limitações:

- depende fortemente do tamanho da população;
- não permite comparar bem territórios com populações muito diferentes;
- pode aumentar apenas porque a população é maior ou mais envelhecida.

### Mortalidade Bruta

A mortalidade bruta é o número de óbitos dividido pela população, geralmente apresentado por 100.000 habitantes.

Interpretação simples:

```text
mortalidade bruta = óbitos / população * 100.000
```

Vantagens:

- fácil de compreender;
- mostra o risco observado na população real;
- útil para descrever carga de mortalidade num território.

Limitações:

- é afectada pela estrutura etária;
- locais mais envelhecidos tendem a ter taxas brutas mais altas;
- não é a melhor métrica para comparar regiões com idades muito diferentes.

### Mortalidade Padronizada ou Ajustada

A mortalidade padronizada, também chamada mortalidade ajustada por idade, tenta remover parte do efeito da estrutura etária. A aplicação usa padronização directa com a População Padrão Europeia de 2013.

Ideia central:

- calcula taxas específicas por idade;
- aplica essas taxas a uma população padrão;
- devolve uma taxa comparável entre locais ou anos com estruturas etárias diferentes.

Vantagens:

- é mais adequada para comparar regiões;
- ajuda a separar diferenças reais de mortalidade de diferenças causadas pelo envelhecimento populacional;
- é preferível para séries temporais quando a estrutura etária muda ao longo do tempo.

Limitações:

- é menos intuitiva do que a mortalidade bruta;
- pode ser instável quando há poucos óbitos em algumas idades;
- não elimina todos os factores de confundimento;
- depende da população padrão escolhida.

Quando usar:

- para comparar `Portugal`, `Norte` e uma região local;
- para comparar uma mesma região ao longo de muitos anos;
- para previsões em que o interesse principal é o padrão de mortalidade ajustado à idade.

### Mortalidade Proporcional

A mortalidade proporcional mostra que percentagem dos óbitos totais pertence a uma causa específica.

```text
mortalidade proporcional = óbitos por causa / óbitos por todas as causas * 100
```

Vantagens:

- mostra o peso relativo de uma causa dentro da mortalidade total;
- é útil para comparar prioridades relativas;
- pode ser informativa quando não se quer trabalhar directamente com população denominadora.

Limitações importantes:

- não mede risco de morrer dessa causa na população;
- pode aumentar porque outras causas diminuíram;
- pode diminuir mesmo que os óbitos dessa causa se mantenham estáveis, se outras causas aumentarem;
- exige dados de `Todas as causas de morte` para o mesmo local, ano e sexo.

Quando a aplicação calcula `Mortalidade Proporcional`, carrega `Todas as causas de morte` como denominador, mesmo que essa causa não tenha sido seleccionada directamente.

### AVPP

`AVPP` significa Anos de Vida Potencialmente Perdidos. A aplicação usa 70 anos como ponto de corte.

A ideia é dar mais peso às mortes que ocorrem em idades mais jovens. Uma morte aos 40 anos contribui mais AVPP do que uma morte aos 68 anos; uma morte depois dos 70 anos não contribui para este indicador.

Vantagens:

- destaca mortalidade prematura;
- ajuda a identificar causas com impacto em idades mais jovens;
- pode complementar taxas de mortalidade, que muitas vezes são dominadas por idades avançadas.

Limitações:

- é uma aproximação porque os dados estão agrupados por bandas etárias;
- depende do ponto de corte escolhido;
- não deve ser comparado como se fosse uma taxa padronizada, salvo se houver uma metodologia adicional para isso.

## 5. Incerteza e Intervalos de Confiança

Quando possível, a aplicação apresenta intervalos de confiança a 95%.

Para mortalidade bruta, usa intervalos baseados em contagens de óbitos. Para mortalidade proporcional, usa uma aproximação binomial. Para AVPP, usa uma aproximação baseada na variação esperada das contagens por idade. Para mortalidade padronizada, usa a rotina de padronização directa disponível no pacote utilizado.

Como interpretar:

- intervalos mais estreitos sugerem estimativas mais precisas;
- intervalos largos são comuns em áreas pequenas, causas raras ou poucos anos;
- intervalos que se sobrepõem não provam que não há diferença, mas aconselham cautela;
- intervalos que não se sobrepõem também não substituem uma análise estatística formal.

Os intervalos ajudam a lembrar que uma taxa observada é uma estimativa, não uma verdade fixa.

## 6. Separador Mortalidade Observada

Use este separador para ver a evolução histórica de uma causa de morte.

Fluxo típico:

1. escolha `Local de residência`;
2. escolha `Causa de Morte`;
3. escolha `Sexo`;
4. escolha `População`;
5. escolha `Taxa`;
6. escolha `Fonte de dados`;
7. escolha os anos a importar;
8. carregue os dados.

Resultados principais:

- gráfico temporal;
- tabela anual;
- resumo da selecção;
- indicação dos indicadores usados como fonte;
- botões para exportar tabelas e imagens.

Quando usar:

- para observar tendências históricas;
- para ver diferenças entre taxa bruta e taxa padronizada;
- para preparar uma série que depois será analisada nos separadores de previsão;
- para verificar se há anos ausentes ou instáveis antes de avançar.

Cuidados:

- confirme se a fonte usada é RDS ou INE em directo;
- em causas raras, pequenas flutuações podem produzir grandes variações nas taxas;
- em áreas pequenas, observe sempre os intervalos de confiança.

## 7. Separador Previsão Guiada

Este separador é para obter uma previsão rápida com menos decisões técnicas. Foi pensado para utilizadores que querem uma projecção exploratória sem configurar manualmente todos os modelos.

Controlos principais:

- `Local de residência`;
- `Nome da Selecção (opcional)`;
- `Causa de Morte`;
- `Sexo`;
- `População`;
- `Taxa`;
- `Fonte de dados`;
- `Janela de ajuste`;
- horizonte da previsão;
- modo de previsão;
- `Como escolher o modelo recomendado` (método de validação);
- `Tamanho do teste (% dos anos)`, quando aplicável.

A `Janela de ajuste` define os anos usados para treinar o modelo. Se usar uma janela mais curta, a previsão fica mais focada na tendência recente. Se usar uma janela mais longa, fica mais influenciada pela história completa.

Modos típicos:

- previsão recomendada;
- comparação entre modelos disponíveis.

### Como é Escolhido o Modelo Recomendado

Por predefinição, o modelo recomendado é escolhido pela sua precisão **fora da amostra**, e não apenas pela qualidade do ajuste à série completa. Assim evita-se favorecer modelos que se ajustam muito bem ao passado mas que preveem mal, um risco real em séries anuais curtas.

Os anos mais recentes da série formam o período de teste. O controlo `Tamanho do teste (% dos anos)` define que percentagem dos anos é usada para esse teste. Há dois esquemas, mais um modo de referência:

- `Validação móvel (recomendada)`: para cada ano do período de teste, o modelo é reajustado com os anos anteriores e avaliado numa previsão a um passo; os erros são combinados. Aproveita melhor as séries curtas e é menos sensível a um único corte.
- `Divisão única treino/teste`: o modelo é ajustado uma vez nos anos iniciais e avaliado no período de teste completo de uma só vez. É mais simples, mas mais sensível ao período de teste escolhido.
- `Ajuste dentro da amostra`: usa apenas a precisão no ajuste à série completa. É o comportamento mais simples e serve sobretudo como referência.

Se a série for demasiado curta para reservar pelo menos três anos de treino e um ano de teste, a aplicação recorre automaticamente ao ajuste dentro da amostra e indica-o no painel de fiabilidade. A validação avalia a previsão a um passo, pelo que reflecte sobretudo o desempenho de curto prazo; horizontes longos continuam a exigir cautela.

Resultados principais:

- gráfico da série observada e valores previstos;
- tabela com anos observados e previstos;
- resumo do modelo recomendado;
- aviso de fiabilidade quando a série é curta, incompleta ou instável;
- botões para exportar tabela e gráfico.

Quando usar:

- para uma primeira projecção;
- para apresentar uma previsão simples;
- para comparar rapidamente se vários modelos dão mensagens semelhantes;
- para detectar se os dados são insuficientes antes de ir para a previsão avançada.

Cuidados:

- uma previsão não é uma meta nem uma previsão oficial;
- horizontes longos, como até 2050, devem ser lidos com cautela;
- se a aplicação mostrar `Erro detectado na previsão`, não interprete os valores como resultados válidos;
- se houver dados em falta, a previsão pode ser bloqueada ou acompanhada por aviso.

## 8. Separador Previsão Avançada

Este separador é para análises mais técnicas. Permite definir a especificação do modelo, comparar métodos, ver diagnósticos, fazer backtesting e explorar quebras estruturais.

Use quando precisa de:

- controlar modelos específicos;
- comparar métodos de previsão;
- avaliar resíduos;
- testar desempenho em períodos de validação;
- procurar mudanças estruturais na série;
- justificar melhor a escolha do modelo.

Áreas principais:

- especificação dos dados e modelos;
- resultados da previsão;
- diagnósticos;
- comparação e métricas de erro;
- quebras estruturais.

### Modelos de Previsão

`ARIMA` modela tendência, diferenças e autocorrelação. É flexível e muitas vezes eficaz em séries temporais, mas pode ser difícil de explicar e pode sobreajustar séries curtas.

`ETS` usa suavização exponencial em modelos de erro, tendência e nível. Costuma funcionar bem em séries suaves, mas pode ter dificuldade com alterações súbitas.

`Random walk with drift` assume que a série continua a tendência média recente. É transparente, mas pode prolongar uma tendência que já não faz sentido epidemiológico.

`Naive` assume que o futuro é igual ao último valor observado. É uma referência simples. Se modelos complexos não forem claramente melhores do que o naive, isso é um sinal de cautela.

`Theta` é um método simples que pode funcionar bem em séries com tendência. A interpretação é menos directa do que a de uma linha de tendência simples.

`TBATS` é um modelo flexível para tendências e sazonalidade complexa. Em mortalidade anual pode ser excessivo, porque não há uma sazonalidade intra-anual explícita na série anual.

`Holt` é uma suavização exponencial com tendência. Pode ser útil quando a série tem uma tendência aproximadamente estável.

`Holt amortecido` reduz a extrapolação indefinida da tendência. Pode ser mais prudente quando uma subida ou descida não deve continuar para sempre ao mesmo ritmo.

### Transformações

A transformação logarítmica pode estabilizar séries positivas e reduzir a influência de valores extremos. A previsão é ajustada na escala transformada e depois reconvertida.

Vantagens:

- pode melhorar a estabilidade de séries com variação proporcional;
- reduz a probabilidade de previsões negativas em taxas positivas.

Limitações:

- torna a interpretação menos directa;
- não resolve problemas de dados escassos;
- pode distorcer séries com muitos valores próximos de zero.

### Diagnósticos

Os diagnósticos ajudam a perceber se o modelo deixou padrões por explicar.

Verifique:

- resíduos ao longo do tempo;
- distribuição dos resíduos;
- autocorrelação dos resíduos;
- avisos de estimação;
- comparação entre observado e ajustado.

Sinais de alerta:

- resíduos com tendência clara;
- autocorrelação forte;
- grandes valores extremos;
- modelos diferentes com previsões muito divergentes;
- série curta ou com muitos anos em falta.

### Backtesting

O backtesting reserva os anos mais recentes como período de teste, treina o modelo nos anos anteriores e compara as previsões com os valores que realmente ocorreram. O controlo `Tamanho do teste (% dos anos)` define que percentagem dos anos entra no teste. A `Abordagem de validação` tem três opções:

- `Métricas do ajuste actual`: usa a precisão no ajuste à série completa (dentro da amostra), sem reservar anos.
- `Divisão única (últimos %)`: ajusta uma vez nos anos iniciais e avalia todo o período de teste de uma só vez.
- `Validação móvel (últimos %)`: reajusta o modelo em cada origem do período de teste e avalia previsões a um passo, combinando os erros. É a predefinição e a mais robusta em séries curtas.

A abordagem escolhida aqui é a mesma que determina o **modelo recomendado** usado por predefinição no separador de resultados e nos diagnósticos, pelo que mudar de abordagem ou de percentagem pode alterar o modelo destacado. Se a série for demasiado curta para reservar treino e teste, a aplicação recorre ao ajuste dentro da amostra.

Vantagens:

- dá uma noção prática de desempenho fora da amostra;
- ajuda a comparar modelos;
- pode revelar modelos que parecem bons no ajuste mas falham na previsão;
- a validação móvel aproveita mais os poucos anos disponíveis do que um único corte.

Limitações:

- há poucos anos disponíveis em muitas séries;
- um único período de teste pode não representar o futuro;
- a validação a um passo reflecte sobretudo o desempenho de curto prazo;
- mudanças excepcionais, como epidemias ou alterações de codificação, podem distorcer o resultado.

### Análise de Quebras

A análise de quebras procura mudanças na estrutura da série, por exemplo alterações no nível médio ou na tendência.

Pode ser útil para levantar hipóteses sobre:

- alterações de codificação ou classificação;
- mudanças epidemiológicas reais;
- impacto de eventos excepcionais;
- transições demográficas;
- mudanças nos sistemas de registo ou cobertura.

Limitações fundamentais:

- uma quebra estatística não prova a causa da mudança;
- séries curtas tornam a detecção menos robusta;
- causas raras podem mostrar quebras por instabilidade aleatória;
- uma quebra deve ser interpretada com conhecimento epidemiológico e histórico.

Se houver uma quebra importante, compare previsões com janelas de ajuste diferentes. Uma previsão usando toda a série pode ser menos adequada do que uma previsão usando apenas o período posterior à quebra.

## 9. Métricas de Erro

As métricas de erro ajudam a comparar modelos. Nenhuma métrica escolhe o modelo perfeito sozinha.

`ME` é o erro médio. Mostra viés médio. Valores positivos ou negativos indicam se o modelo tende a subestimar ou sobrestimar.

Vantagem: mostra direcção do erro.

Limitação: erros positivos e negativos podem anular-se.

`RMSE` é a raiz do erro quadrático médio. Penaliza mais erros grandes.

Vantagem: útil quando erros grandes são especialmente importantes.

Limitação: pode ser dominado por poucos anos extremos.

`MAE` é o erro absoluto médio. Mede o erro médio em unidades da série.

Vantagem: mais robusto a extremos do que RMSE.

Limitação: não penaliza erros grandes tão fortemente.

`MAPE` é o erro percentual absoluto médio.

Vantagem: intuitivo por estar em percentagem.

Limitação: pode ser problemático quando os valores observados são pequenos ou próximos de zero.

`MASE` é o erro absoluto médio escalado. Compara o erro do modelo com um modelo naive.

Vantagem: permite comparar desempenho entre séries com escalas diferentes.

Limitação: é menos intuitivo para utilizadores não técnicos.

Sugestão prática:

- veja RMSE e MAE em conjunto;
- use MAPE com cautela em causas raras;
- confirme se o modelo escolhido também faz sentido no gráfico;
- não aceite um modelo só porque uma métrica ficou ligeiramente melhor.

## 10. Separador Métricas Anuais

Este separador compara métricas num único ano. Mostra sempre:

- `Portugal`;
- `Norte`;
- a localização adicional seleccionada.

Pode seleccionar uma ou várias causas de morte. A tabela é ordenada do maior para o menor valor segundo a localização adicional. Isto ajuda a perceber quais as causas com maior peso local, mantendo Portugal e Norte como comparadores.

Métricas disponíveis:

- `Óbitos`;
- `Mortalidade Bruta`;
- `Mortalidade Padronizada`;
- `Mortalidade Proporcional`;
- `AVPP`.

Quando usar:

- para comparar prioridades num ano específico;
- para ver se o perfil local difere de Portugal ou Norte;
- para seleccionar causas que merecem análise temporal ou previsão;
- para preparar quadros de apresentação.

Cuidados:

- para `Mortalidade Proporcional`, lembre-se de que o denominador é `Todas as causas de morte`;
- causas raras podem surgir com taxas instáveis;
- ordenação por valor local não significa importância causal ou evitabilidade.

## 11. Separador Disponibilidade de Dados

Este separador mostra a cobertura dos ficheiros RDS. É útil antes de carregar análises pesadas.

Estados possíveis:

- `Disponível`: todos os anos, causas e áreas pedidas estão presentes;
- `Parcial`: há dados para parte da selecção, mas não para tudo;
- `Indisponível`: não há dados suficientes nos RDS para essa selecção.

Quando usar:

- antes de seleccionar muitas causas;
- antes de usar anos antigos do indicador `0008206`;
- quando uma análise devolve aviso de dados incompletos;
- para decidir se vale a pena usar `INE em directo`.

Esta verificação é feita ao nível do inventário. A análise final ainda pode falhar se os ficheiros existirem mas não tiverem as linhas esperadas.

## 12. Exportação de Resultados

Os separadores têm botões para exportar:

- tabelas em CSV;
- gráficos em PNG.

Antes de exportar:

- confirme a fonte de dados;
- confirme ano, local, causa, sexo e população;
- verifique se há avisos;
- evite exportar previsões com erro detectado ou dados incompletos sem mencionar essa limitação.

## 13. Fluxos de Trabalho Recomendados

### Explorar uma Tendência Observada

1. Abra `Mortalidade Observada`.
2. Escolha `Ficheiros RDS`.
3. Seleccione local, causa, sexo, população e taxa.
4. Use uma janela de anos ampla.
5. Observe gráfico, tabela e intervalos.
6. Exporte se a série estiver coerente.

### Fazer uma Previsão Simples

1. Abra `Previsão Guiada`.
2. Escolha a mesma combinação que pretende analisar.
3. Defina a janela de ajuste.
4. Escolha horizonte.
5. Veja a previsão recomendada e os avisos.
6. Compare modelos se houver dúvida.

### Fazer uma Previsão Técnica

1. Abra `Previsão Avançada`.
2. Defina cuidadosamente fonte, local, causa, sexo, população e taxa.
3. Escolha modelos candidatos.
4. Analise métricas de erro.
5. Veja diagnósticos.
6. Faça backtesting se houver anos suficientes.
7. Veja análise de quebras antes de aceitar uma janela de ajuste longa.

### Comparar Causas num Ano

1. Abra `Métricas Anuais`.
2. Escolha o ano.
3. Seleccione várias causas.
4. Escolha a métrica.
5. Compare Portugal, Norte e a localização adicional.
6. Use a ordenação local para identificar causas prioritárias para análise posterior.

### Resolver Problemas de Dados

1. Abra `Disponibilidade de Dados`.
2. Escolha anos, locais e causas.
3. Confirme se os RDS cobrem a selecção.
4. Se estiver `Parcial` ou `Indisponível`, reduza a selecção ou use `INE em directo`.

## 14. Boas Práticas de Interpretação

- Compare mortalidade bruta e padronizada quando a estrutura etária for relevante.
- Use mortalidade padronizada para comparações geográficas sempre que possível.
- Não interprete mortalidade proporcional como risco populacional.
- Em áreas pequenas, privilegie padrões consistentes em vez de um único ano.
- Em causas raras, observe sempre os intervalos e a estabilidade temporal.
- Em previsões longas, destaque que são extrapolações estatísticas.
- Se houver quebras estruturais, faça análises de sensibilidade com janelas de ajuste diferentes.
- Documente sempre a fonte de dados usada: RDS ou INE em directo.
- Verifique se os dados podem ter sido revistos pelo INE.

## 15. Limitações Gerais

A aplicação não ajusta para factores individuais, comorbilidades, privação socioeconómica, exposição ambiental, acesso a cuidados ou alterações de diagnóstico.

As taxas padronizadas ajudam a comparar estruturas etárias diferentes, mas não resolvem todos os problemas de comparabilidade.

As previsões assumem que padrões históricos carregam informação sobre o futuro. Essa suposição pode falhar quando há alterações epidemiológicas, tecnológicas, sociais, ambientais ou de codificação.

Use a aplicação como apoio à análise, não como resposta final.
