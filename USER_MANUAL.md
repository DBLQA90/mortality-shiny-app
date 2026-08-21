# Manual do Utilizador

Este manual explica como usar a aplicação e como interpretar os principais conceitos epidemiológicos e estatísticos apresentados nos separadores. A aplicação é uma ferramenta exploratória não oficial: ajuda a analisar dados de mortalidade do INE, mas não substitui validação epidemiológica, análise clínica, revisão metodológica ou produtos estatísticos oficiais.

Para detalhes técnicos das fórmulas, indicadores e implementação, ver [METHODOLOGY.md](METHODOLOGY.md).

## 1. Visão Geral

A aplicação permite:

- consultar mortalidade observada por ano, local, causa de morte, sexo e grupo populacional;
- comparar Portugal, Norte e uma localização adicional num ano específico;
- calcular nove métricas: óbitos, mortalidade bruta, mortalidade padronizada (directa),
  SMR e taxa padronizada indirecta, mortalidade proporcional, AVPP, óbitos infantis e
  mortalidade infantil;
- agregar três ou cinco anos numa só janela, para estabilizar concelhos pequenos e causas raras;
- escolher a definição das regiões, NUTS 2013 ou NUTS 2024;
- gerar previsões exploratórias até 30 anos no futuro;
- comparar modelos de previsão;
- avaliar diagnósticos, erros de previsão e possíveis quebras estruturais;
- verificar a disponibilidade de dados nos ficheiros RDS antes de carregar uma análise;
- exportar tabelas em CSV e gráficos em PNG.

A ordem dos separadores é:

1. `Introdução`
2. `Mortalidade Observada`
3. `Previsão Guiada`
4. `Previsão Avançada`
5. `Métricas Anuais`
6. `Disponibilidade de Dados`
7. `Glossário`

O separador `Introdução` é a página inicial: explica em linguagem simples o que é uma previsão, como começar em três passos e qual separador usar. É o ponto de partida recomendado para quem abre a aplicação pela primeira vez.

Acima dos separadores há um controlo que se aplica a toda a aplicação, `Definição das regiões`. Escolhe se as regiões seguem a definição NUTS 2013 ou NUTS 2024. Não altera os dados lidos, apenas a forma como os municípios são agrupados — ver a secção *Definição das Regiões (NUTS)*.

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

Define a área geográfica usada na análise. A lista tem três tipos de entrada:

- `Portugal`, o total nacional;
- as regiões: `Continente`, as duas regiões autónomas e as regiões NUTS II;
- os 308 municípios.

Pode seleccionar uma única localização ou várias. Quando selecciona mais do que uma, a aplicação soma os óbitos e a população dessas localizações antes de calcular as taxas. Isto é útil para criar uma área agregada, por exemplo uma área de influência, um conjunto de concelhos ou uma ULS.

O campo `Nome da Selecção (opcional)` permite dar um nome a essa agregação. Se for deixado em branco, a aplicação usa uma designação automática.

**As regiões são sempre somadas a partir dos seus municípios.** Quando escolhe `Alentejo`, a aplicação não lê a linha «Alentejo» publicada pelo INE: soma os municípios que hoje pertencem ao Alentejo, e aplica essa mesma lista a todos os anos da série. Isto tem uma razão e uma consequência.

A razão é que o INE mudou as fronteiras das regiões em 2024, e a Lezíria do Tejo passou do Alentejo para a nova região Oeste e Vale do Tejo. Ler as linhas regionais do INE ao longo da série significaria comparar dois Alentejos diferentes: os óbitos de 2022 seriam 11.327 numa definição e 7.898 na outra. Somando sempre os mesmos municípios, a série mantém-se contínua e comparável.

A consequência é que os totais regionais **não coincidem exactamente com os publicados pelo INE**. A diferença é pequena — vem dos óbitos que o INE não conseguiu atribuir a nenhum município — mas existe. Quem precisar de reproduzir um número publicado deve usar a fonte oficial, não esta aplicação.

Duas notas práticas:

- Seleccionar uma região **e** um município que lhe pertence não duplica nada. Como a região é expandida numa lista de municípios sem repetições, `Alentejo` + `Beja` dá o Alentejo. A aplicação avisa quando isso acontece, porque provavelmente não era o que pretendia.
- Seleccionar `Portugal` com qualquer outra coisa **duplica**, porque `Portugal` é lido como uma linha própria e não expandido. Também aqui há aviso.

### Definição das Regiões (NUTS)

O controlo `Definição das regiões`, no topo de qualquer página, escolhe qual das duas versões da nomenclatura NUTS agrupa os municípios:

| | Regiões | Lisboa | Lezíria do Tejo, Oeste, Médio Tejo |
|---|---:|---|---|
| **NUTS 2013** | 7 | uma região, `Área Metropolitana de Lisboa` | dentro de `Alentejo` e `Centro` |
| **NUTS 2024** | 9 | `Grande Lisboa` + `Península de Setúbal` | formam `Oeste e Vale do Tejo` |

As duas definições cobrem exactamente os mesmos 308 municípios. Mudar de definição **não altera os dados lidos**, apenas a forma como são agrupados, e qualquer das duas dá uma série contínua em todos os anos.

O ponto a que deve prestar atenção: **seis nomes existem nas duas definições e significam coisas diferentes em cada uma**. O `Centro` de 2013 tem 100 municípios; o de 2024 tem 77. O `Alentejo` de 2013 inclui a Lezíria do Tejo; o de 2024 não. É por isso que o controlo está no topo da página e não escondido num separador: a definição activa tem de estar à vista sempre que lê um número regional.

Se mudar de definição com uma região seleccionada que não existe na outra — `Grande Lisboa` não existe em NUTS 2013 — essa região é retirada da selecção e a aplicação diz quais foram retiradas. Os municípios seleccionados nunca são afectados.

Quando usar NUTS 2013:

- para acompanhar documentos, relatórios ou séries publicadas antes da revisão de 2024;
- quando precisa da Área Metropolitana de Lisboa como uma única região;
- quando quer reproduzir aproximadamente uma linha regional publicada pelo INE: nesta definição, e nos anos em que o arquivo tem a linha do INE para comparação, as somas municipais coincidem com ela.

Quando usar NUTS 2024 (predefinição):

- para trabalho corrente e para qualquer coisa que vá ser comparada com publicações recentes;
- quando quer distinguir a Grande Lisboa da Península de Setúbal.

### Continente, Açores e Madeira

`Continente`, `Região Autónoma dos Açores` e `Região Autónoma da Madeira` são o nível NUTS I, acima das regiões. Estão disponíveis nas duas definições e comportam-se como qualquer outra região: são somados a partir dos seus municípios — 278 no Continente, 19 nos Açores, 11 na Madeira.

Os três repartem o país sem sobreposição nem falha. Para 2021, 119.589 + 2.366 + 2.875 = 124.830 óbitos, que é exactamente o total nacional publicado pelo INE.

Seleccionar `Continente` junto com uma região que lhe pertence, como `Norte`, não duplica: o Continente já inclui o Norte, e o resultado é apenas o Continente. A aplicação avisa.

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

### Gráficos Interactivos

Os principais gráficos de taxas e de previsão são interactivos. Passe o rato sobre um ponto para ver o ano e o valor exactos. Pode também aproximar (zoom) uma zona do gráfico arrastando o rato e voltar ao início com os botões da barra de ferramentas que aparece no canto do gráfico. Para guardar uma imagem fixa, use o botão `Descarregar gráfico (PNG)`.

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

### SMR (Padronização Indirecta)

O `SMR` compara os óbitos que um local teve com os que teria tido se tivesse as taxas por idade de uma referência. A referência vale 100.

```text
SMR = óbitos observados / óbitos esperados * 100
```

Os óbitos esperados calculam-se aplicando as taxas por idade da referência à estrutura etária do local. Se o Alentejo tem muitas pessoas idosas, espera-se que tenha muitos óbitos mesmo que a mortalidade em cada idade seja igual à do país; o SMR pergunta se teve mais ou menos do que isso.

Leitura:

- `SMR = 100`: o local tem a mortalidade que se esperaria dada a sua estrutura etária;
- `SMR = 120`: 20% mais óbitos do que o esperado;
- `SMR = 85`: 15% menos óbitos do que o esperado.

**A diferença face à padronização directa** está em que estrutura é usada. A padronização directa aplica as taxas *do local* a uma população padrão externa (a europeia de 2013): precisa de uma taxa estimável em cada banda etária, e num concelho pequeno muitas bandas têm zero óbitos, o que torna a taxa instável ou impossível de calcular. A padronização indirecta faz o inverso — aplica as taxas *da referência*, que são estáveis porque vêm de uma população grande, à estrutura do local. Só precisa do total de óbitos observados no local, não de uma taxa por cada idade.

Por isso o SMR é a métrica indicada para concelhos pequenos e causas raras, onde a padronizada directa é instável ou nem sequer estimável.

Vantagens:

- estável em populações pequenas, onde a padronização directa falha;
- leitura imediata: 100 é a referência;
- usa toda a informação disponível sem exigir taxas por idade no local.

Limitações:

- **dois SMR não se comparam bem entre si.** Cada um é calculado contra a estrutura etária do seu próprio local, por isso «Beja 120» e «Braga 110» não significam que Beja tenha 9% mais mortalidade do que Braga. Cada SMR compara-se com a referência, não com outro SMR;
- depende da referência escolhida;
- como qualquer indicador local, tem intervalos largos quando há poucos óbitos.

`Referência da padronização indirecta` escolhe contra quem a comparação é feita. `Portugal` é a predefinição, e faz com que 100 signifique sempre «igual à média nacional».

### Taxa Padronizada Indirecta

A mesma conta do SMR, expressa como taxa por 100.000 em vez de como índice. Obtém-se multiplicando o SMR pela taxa bruta da referência.

Serve para quem prefere ler uma taxa a ler um índice. As vantagens, limitações e cuidados são exactamente os do SMR — em particular, duas taxas padronizadas indirectas de locais diferentes continuam a não ser directamente comparáveis entre si.

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

Uma nota sobre a banda `0 - 4 anos`. Como a idade de cada morte é aproximada pelo ponto médio da sua banda, todas as mortes entre os 0 e os 4 anos contariam 67,5 anos perdidos. Mas a maioria delas são mortes no primeiro ano de vida — em 2024, 254 dos 286 óbitos nacionais nessa banda — e essas perdem quase os 70 anos completos. A aplicação separa a banda em `< 1 ano` e `1 - 4 anos` usando as contagens de óbitos infantis, o que corrige a subestimação. O efeito é pequeno no total (+0,17% em todas as causas) e maior nas causas perinatais e congénitas (+2,9% e +1,7%).

### Óbitos Infantis (< 1 ano)

Número de mortes antes do primeiro ano de vida, sem denominador.

É a métrica a usar à escala municipal. A taxa de mortalidade infantil (abaixo) precisa de um denominador que, num concelho pequeno, é minúsculo; a contagem diz o que aconteceu e não pode ser mal lida como se fosse comparável entre locais.

Duas particularidades:

- ao contrário das outras contagens, **não é convertida em média anual** quando agrega vários anos. Dois óbitos em três anos apareceriam como «1», uma fracção arredondada, quando o que aconteceu foram dois óbitos. Fica o total do período, que é também exactamente o numerador da taxa apresentada ao lado;
- cobre 1991-2025, mais do que a taxa, porque só precisa do numerador.

### Mortalidade Infantil (por 1.000 nados-vivos)

Óbitos com menos de 1 ano por 1.000 nados-vivos.

```text
mortalidade infantil = óbitos com menos de 1 ano / nados-vivos * 1.000
```

**O denominador são os nados-vivos, não a população.** Esta é a única taxa da aplicação cujo denominador não é uma população, e a razão é simples: nenhum indicador de população do INE tem uma banda etária «menos de 1 ano», por isso «óbitos infantis por população com menos de 1 ano» não é calculável de todo. Os nados-vivos são a convenção internacional, e correspondem melhor ao grupo em risco: quem pode morrer no primeiro ano de vida é quem nasceu.

Também por isso a escala é por 1.000 e não por 100.000 como as restantes taxas.

Cobertura: 1995-2025. Em 2025 o INE publica os óbitos com menos de 1 ano apenas no total, sem separar por causa nem por sexo, e a aplicação recusa pedidos mais detalhados nesse ano em vez de responder zero.

Vantagens:

- é o indicador clássico de saúde materno-infantil, comparável internacionalmente;
- reconcilia com a série publicada pelo INE em todos os anos.

Limitações:

- **à escala municipal é extremamente esparso.** Barrancos teve 9 nados-vivos em 2024. Sem óbitos infantis, a taxa é 0,0 — mas o limite superior do intervalo é 409,9 por 1.000. Um único óbito teria dado mais de 100;
- por isso os valores calculados sobre menos de 1.000 nados-vivos são assinalados com `*` (ver a secção seguinte);
- agregar vários anos ajuda, mas não cria acontecimentos que não houve.

### Que Métrica Escolher

A aplicação tem nove métricas, e a escolha certa depende mais da pergunta e do tamanho do local do que de preferência pessoal.

| A sua pergunta | Métrica |
|---|---|
| Quantas pessoas morreram? | `Óbitos` |
| Qual é o risco observado nesta população? | `Mortalidade Bruta` |
| Como se compara este local com outro, sendo as idades diferentes? | `Mortalidade Padronizada (directa)` |
| E se o local for pequeno, com poucos óbitos? | `SMR` |
| Que peso tem esta causa no total de mortes? | `Mortalidade Proporcional` |
| Que causas matam mais cedo? | `AVPP` |
| Como está a saúde materno-infantil? | `Mortalidade Infantil` (país, região) ou `Óbitos infantis` (concelho) |

Regras práticas:

- **Nunca compare territórios com mortalidade bruta** se as estruturas etárias diferirem, e em Portugal diferem quase sempre. Um concelho do interior tem taxa bruta alta sobretudo porque é envelhecido.
- **Padronizada directa para regiões, SMR para concelhos.** A fronteira não é rígida, mas se o intervalo de confiança da padronizada directa for muito largo, ou se o valor não for calculável, mudou de regime e deve usar o SMR.
- **Compare cada SMR com 100, não com outro SMR.**
- **Mortalidade proporcional não é risco.** Uma causa pode subir de peso apenas porque outras desceram.
- **Se a série tiver muito ruído, agregue 3 ou 5 anos** antes de mudar de métrica. Muitas vezes é o tamanho do denominador, não a métrica, que está a causar o problema.

## 5. Incerteza e Intervalos de Confiança

Quando possível, a aplicação apresenta intervalos de confiança a 95%.

Para mortalidade bruta, usa intervalos baseados em contagens de óbitos. Para mortalidade proporcional, usa uma aproximação binomial. Para AVPP, usa uma aproximação baseada na variação esperada das contagens por idade. Para mortalidade padronizada, usa a rotina de padronização directa disponível no pacote utilizado. Para `SMR` e taxa padronizada indirecta, usa o método de Byar. Para mortalidade infantil, um intervalo de Poisson exacto sobre a contagem de óbitos, escalado pelos nados-vivos.

Como interpretar:

- intervalos mais estreitos sugerem estimativas mais precisas;
- intervalos largos são comuns em áreas pequenas, causas raras ou poucos anos;
- intervalos que se sobrepõem não provam que não há diferença, mas aconselham cautela;
- intervalos que não se sobrepõem também não substituem uma análise estatística formal.

Os intervalos ajudam a lembrar que uma taxa observada é uma estimativa, não uma verdade fixa.

### Um intervalo largo não é um defeito

É a leitura errada mais frequente. Quando um concelho pequeno apresenta um intervalo enorme, a aplicação não está a falhar: está a dizer com honestidade que, com aquele número de acontecimentos, não é possível saber mais. Barrancos, 2024, mortalidade infantil: `0 (0; 409,88)`. O zero é real — não houve óbitos infantis. O 409,88 também é real — com 9 nados-vivos, um único óbito daria mais de 100 por 1.000.

O que fazer perante um intervalo largo:

- **não conclua que o local é melhor ou pior** do que outro se os intervalos se sobrepõem largamente;
- **agregue 3 ou 5 anos**, que é o instrumento que a aplicação oferece exactamente para isto;
- **mude para o SMR** se o problema for a instabilidade da padronização directa;
- **use a contagem** em vez da taxa quando o denominador é minúsculo;
- **não elimine o valor** do relatório por ser incerto: apresente-o com o intervalo.

### O asterisco na mortalidade infantil

Um valor de mortalidade infantil marcado com `*` foi calculado sobre **menos de 1.000 nados-vivos no período**. A marca aparece na tabela, no ficheiro CSV exportado e no gráfico, sempre acompanhada de uma nota de rodapé.

O limiar não é um teste de significância, é uma afirmação sobre resolução: abaixo de 1.000 nados-vivos, um único óbito adicional desloca a taxa em mais de uma unidade inteira por 1.000 — mais do que toda a taxa nacional, que ronda 3. Ordenar concelhos por esse valor, ou ler uma variação entre anos, é ler ruído.

A maioria dos concelhos portugueses fica abaixo do limiar, e isso é precisamente o ponto: a marca descreve o caso normal, não um punhado de excepções. **Nada é escondido** — o valor é apresentado, é exacto, e o intervalo já indica a incerteza. O asterisco existe apenas para impedir que alguém, ao percorrer a tabela, trate o número como comparável.

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

### O que é uma previsão

Uma previsão, ou projecção, é uma estimativa de como uma taxa poderá evoluir no futuro, a partir do padrão dos anos anteriores. É importante ter presente que:

- não é uma certeza: é um cenário possível, não o que vai necessariamente acontecer;
- não é uma meta nem um número oficial;
- a incerteza aumenta com o tempo, pelo que os primeiros anos são mais fiáveis do que, por exemplo, daqui a 20 ou 30 anos.

A aplicação aprende o padrão a partir dos anos observados (a janela de ajuste) e prolonga-o para o futuro, apresentando também uma banda de incerteza. Interprete os resultados como apoio à exploração, e não como conclusões definitivas. Os termos usados aqui estão explicados no separador `Glossário`.

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
- em `Mostrar opções avançadas` (opcional): o método de validação (`Como escolher o modelo recomendado`) e o `Tamanho do teste (% dos anos)`.

Por predefinição, a aplicação usa boas escolhas automáticas, pelo que basta escolher local, causa e horizonte e clicar em `Gerar previsão`. As opções de validação só são necessárias para afinar como o modelo é escolhido e ficam escondidas até activar `Mostrar opções avançadas`.

A `Janela de ajuste` define os anos usados para treinar o modelo. Se usar uma janela mais curta, a previsão fica mais focada na tendência recente. Se usar uma janela mais longa, fica mais influenciada pela história completa.

Modos típicos:

- previsão recomendada;
- comparação entre modelos disponíveis.

### Como é Escolhido o Modelo Recomendado

Estas opções ficam em `Mostrar opções avançadas` e não são necessárias para uma previsão simples. Por predefinição, o modelo recomendado é escolhido pela sua precisão **fora da amostra**, e não apenas pela qualidade do ajuste à série completa. Assim evita-se favorecer modelos que se ajustam muito bem ao passado mas que preveem mal, um risco real em séries anuais curtas.

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

A transformação logarítmica soma um pequeno valor de deslocamento (offset) antes de aplicar o logaritmo, para permitir anos com taxa zero. Esse offset corresponde a metade da menor taxa positiva da série de ajuste e é apresentado na etiqueta da transformação (por exemplo, na tabela de especificação do modelo avançado). Quando a série contém zeros, a aplicação assinala que o offset influencia a projecção e os intervalos, porque é nesse caso que uma pequena constante aditiva tem maior efeito. Nessas situações, compare com `Sem transformação` para perceber a sensibilidade do resultado.

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
- `Mortalidade Padronizada (directa, ESP 2013)`;
- `SMR (padronização indirecta, referência = 100)`;
- `Taxa Padronizada Indirecta (por 100.000)`;
- `Mortalidade Proporcional`;
- `AVPP`;
- `Óbitos infantis (< 1 ano)`;
- `Mortalidade Infantil (por 1.000 nados-vivos)`.

### Agregação Plurianual

O controlo `Agregação plurianual` permite calcular a métrica sobre 1, 3 ou 5 anos centrados no ano escolhido. Serve para estabilizar concelhos pequenos e causas raras, onde um único ano tem demasiado poucos acontecimentos para ser lido.

A agregação soma os óbitos **e** a população dos anos incluídos. O denominador passa a ser em *pessoas-ano*: cinco anos de um concelho com 10.000 habitantes são 50.000 pessoas-ano. Como numerador e denominador crescem juntos, a taxa continua a ser por ano e é directamente comparável com um valor não agregado — não é uma soma, é uma média ponderada.

As contagens comportam-se de outra forma. `Óbitos` e `AVPP` são convertidos em **média anual**, porque somar três anos de óbitos triplicaria o número e leria como uma triplicação da mortalidade. `Óbitos infantis` é a excepção deliberada e fica como total do período, pela razão explicada na secção dessa métrica.

Se a janela ultrapassar os anos disponíveis, é truncada, e o período efectivamente usado aparece indicado nos resultados. Uma janela de 5 anos centrada em 2024 usa 2022-2024.

O que ganha e o que perde:

- **ganha** intervalos mais estreitos e séries legíveis em locais pequenos;
- **perde** a capacidade de ver variações anuais reais, que ficam suavizadas.

A agregação aplica-se apenas a este separador. As previsões usam sempre séries anuais não agregadas, porque uma média móvel introduz autocorrelação que invalidaria os intervalos de previsão.

### Referência da Padronização Indirecta

Aparece quando escolhe `SMR` ou `Taxa Padronizada Indirecta`. Define contra quem a comparação é feita: as taxas por idade desta referência são aplicadas à estrutura etária de cada local para calcular os óbitos esperados.

`Portugal` é a predefinição, e faz com que 100 signifique «igual à média nacional».

Quando usar:

- para comparar prioridades num ano específico;
- para ver se o perfil local difere de Portugal ou Norte;
- para seleccionar causas que merecem análise temporal ou previsão;
- para preparar quadros de apresentação.

Cuidados:

- para `Mortalidade Proporcional`, lembre-se de que o denominador é `Todas as causas de morte`;
- causas raras podem surgir com taxas instáveis — considere agregar 3 ou 5 anos, ou mudar para `SMR`;
- ordenação por valor local não significa importância causal ou evitabilidade;
- a definição NUTS activa, no topo da página, muda o que `Centro` ou `Alentejo` significam;
- as taxas exigem denominador: os anos sem população publicada são recusados com uma explicação, e as contagens, a mortalidade proporcional e os AVPP continuam disponíveis.

## 11. Separador Disponibilidade de Dados

Este separador mostra a cobertura dos ficheiros RDS. É útil antes de carregar análises pesadas.

Estados possíveis:

- `Disponível`: todos os anos, causas e áreas pedidas estão presentes;
- `Parcial`: há dados para parte da selecção, mas não para tudo;
- `Indisponível`: não há dados suficientes nos RDS para essa selecção.

O quadro-resumo indica também a cobertura dos quatro conjuntos de dados: população, óbitos por causa, nados-vivos e óbitos com menos de 1 ano. A consulta detalhada, por ano, área e causa, aplica-se aos dois primeiros; os dois conjuntos usados pela mortalidade infantil estão descritos apenas pelo intervalo de anos que cobrem.

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

### Cobertura dos Dados

Os dois lados de uma taxa não terminam no mesmo ano, e é a população que vai à frente:

| Conjunto | Anos | Nota |
|---|---|---|
| População | 1991-2025 | |
| Óbitos por causa | 1991-2024 | limite de tudo excepto a mortalidade infantil |
| Nados-vivos | 1995-2025 | |
| Óbitos com menos de 1 ano | 1980-2025 | 2025 sem detalhe por causa nem por sexo |

Em consequência, **2025 só está disponível para a mortalidade infantil** e é recusado nas restantes métricas com uma mensagem que explica porquê. Existem óbitos de 2025 publicados pelo INE, mas sem a dimensão «causa de morte», por isso não servem esta aplicação.

Todas as métricas funcionam em 2024, incluindo as taxas.

### A Revisão da População de 2021

O INE publica duas estimativas de população que se sobrepõem e **não coincidem**. A série mais recente, em NUTS 2024, revê os valores em alta de forma crescente:

| Portugal | série revista | série anterior | |
|---|---:|---:|---:|
| 2021 | 10.599.117 | 10.421.117 | +1,71% |
| 2022 | 10.929.704 | 10.516.621 | +3,93% |
| 2023 | 11.204.347 | 10.639.726 | +5,31% |

A aplicação usa a série revista em toda a sua extensão, a partir de 2021. A alternativa — usá-la apenas nos anos que faltavam — colocaria um degrau de 5,3% entre 2023 e 2024, que se leria como uma queda real da mortalidade.

Três consequências práticas:

- **todas as taxas a partir de 2021 mudaram.** A mortalidade bruta de Portugal em 2023 é 1.056 por 100.000, e não 1.112 como na base anterior. Se tem números anteriores em circulação, não vão coincidir;
- **existe um degrau em 2020/2021**, de cerca de 1,7%. Uma série de taxas que atravesse esses anos desce ligeiramente por mudança de denominador, não por mudança de mortalidade. A aplicação avisa quando isso acontece;
- **nada antes de 2021 foi revisto.** Verificado contra a série longa nacional do INE, que vai de 1970 a 2025: os valores coincidem ao indivíduo em 1995, 2000, 2005, 2010, 2014, 2016, 2018, 2019 e 2020.

Use a aplicação como apoio à análise, não como resposta final.

## 16. Glossário

Explicações simples dos termos usados na aplicação. O separador `Glossário` apresenta esta mesma lista dentro da aplicação.

### Conceitos de mortalidade

- **Óbitos:** número absoluto de mortes na selecção (ano, local, causa, sexo e idade).
- **Taxa bruta:** número de mortes por 100.000 habitantes. É simples, mas depende muito da idade da população.
- **Taxa padronizada:** taxa ajustada à idade, que permite comparar de forma justa locais com populações mais jovens ou mais envelhecidas. Usa a População Padrão Europeia de 2013.
- **População padrão (ESP 2013):** estrutura etária de referência comum, aplicada na padronização para que as comparações não sejam distorcidas pela idade.
- **Mortalidade proporcional:** percentagem das mortes de uma causa face ao total de mortes, no mesmo ano, sexo e local.
- **AVPP (anos de vida potencialmente perdidos):** medida do impacto da morte prematura; soma os anos que faltavam até aos 70 em cada morte antes dessa idade e dá mais peso às mortes em idades jovens.
- **Mortalidade prematura:** mortes antes de uma certa idade (aqui, antes dos 75 anos), muitas vezes consideradas potencialmente evitáveis.
- **Padronização directa:** aplica as taxas por idade *do local* a uma população padrão externa. Dá a `Taxa padronizada`. Precisa de uma taxa estimável em cada idade, por isso é instável em locais pequenos.
- **Padronização indirecta:** aplica as taxas por idade *da referência* à estrutura etária do local, para calcular quantos óbitos seriam de esperar. Dá o `SMR`. É estável em locais pequenos porque só precisa do total de óbitos observados.
- **SMR:** óbitos observados a dividir pelos esperados, vezes 100. A referência vale 100; 120 são 20% mais óbitos do que o esperado. Compare cada SMR com 100, nunca com outro SMR.
- **Óbitos esperados:** os óbitos que o local teria tido com as taxas por idade da referência e a sua própria estrutura etária.
- **Referência (padronização indirecta):** o território cujas taxas por idade servem de termo de comparação; normalmente `Portugal`.
- **Agregação plurianual:** calcular a métrica sobre 3 ou 5 anos em vez de 1, para estabilizar locais pequenos e causas raras.
- **Pessoas-ano:** o denominador de uma taxa agregada. Cinco anos de um concelho com 10.000 habitantes são 50.000 pessoas-ano, o que mantém a taxa por ano e comparável com um ano isolado.
- **Mortalidade infantil:** óbitos antes do primeiro ano de vida por 1.000 nados-vivos.
- **Nados-vivos:** nascimentos com vida. São o denominador da mortalidade infantil, em vez da população, porque nenhum indicador de população tem uma banda «menos de 1 ano» e porque correspondem melhor ao grupo em risco.
- **Asterisco (`*`):** marca uma taxa de mortalidade infantil calculada sobre menos de 1.000 nados-vivos, em que um único óbito desloca o valor em mais de uma unidade por 1.000. O valor é exacto; a marca avisa que não é comparável.
- **Intervalo de confiança:** margem de incerteza à volta de um valor estimado. Um intervalo de 95% indica uma gama de valores plausíveis; é a zona sombreada nos gráficos. Um intervalo largo não é um defeito: é o que há a dizer quando os acontecimentos são poucos.

### Conceitos de previsão

- **Previsão (projecção):** estimativa de como uma taxa poderá evoluir no futuro, a partir do padrão dos anos anteriores. Não é uma certeza nem uma meta.
- **Horizonte:** quantos anos para o futuro a previsão vai. Quanto maior, maior a incerteza.
- **Janela de ajuste (treino):** os anos usados para o modelo aprender o padrão da série.
- **Teste / validação:** anos recentes reservados para avaliar quão bem o modelo prevê, antes de confiar na projecção futura.
- **Validação móvel:** forma de validação que repete a previsão a partir de várias origens e combina os erros. É a mais fiável em séries curtas.
- **Divisão única (treino/teste):** forma de validação que reserva os últimos anos uma só vez para testar o modelo.
- **Ajuste dentro da amostra:** avaliação usando o ajuste à série completa, sem reservar anos. É menos exigente e serve apenas como referência.
- **Retroteste (backtesting):** testar a previsão contra anos que realmente já aconteceram.
- **Modelo:** método matemático que descreve o padrão da série para o projectar (por exemplo ARIMA, ETS, Holt, Naive). Na Previsão Guiada, a aplicação escolhe um por si.
- **Transformação log:** passo opcional que estabiliza séries positivas e evita previsões negativas; a previsão é feita na escala transformada e depois reconvertida.
- **Métricas de erro (RMSE, MAE, MAPE, MASE):** números que medem quão longe as previsões ficam dos valores reais; servem para comparar modelos. Valores mais baixos são melhores.
- **Quebra estrutural:** mudança no padrão da série (por exemplo no nível ou na tendência), que pode dever-se a alterações reais, de codificação ou de registo.
- **Resíduos e diagnósticos:** os resíduos são as diferenças entre o observado e o ajustado; os diagnósticos (ACF, PACF, Ljung-Box) ajudam a verificar se o modelo captou bem o padrão.

### Dados

- **INE:** Instituto Nacional de Estatística, a fonte oficial dos dados de mortalidade e população.
- **Indicador:** conjunto de dados específico do INE (por exemplo, óbitos por causa), identificado por um código.
- **Ficheiros RDS:** dados já preparados e guardados no repositório, que a aplicação lê rapidamente sem consultar o INE em directo.
- **Fonte de dados:** a escolha entre ler os Ficheiros RDS (rápido) ou consultar o INE em directo (mais lento, para dados não incluídos).

### Geografia

- **NUTS:** a nomenclatura estatística das regiões. A aplicação usa dois níveis: NUTS I (`Continente`, Açores, Madeira) e NUTS II (as regiões).
- **Definição das regiões (NUTS 2013 / NUTS 2024):** as duas versões da nomenclatura que a aplicação oferece, no controlo do topo da página. Agrupam os mesmos 308 municípios de formas diferentes.
- **Área Metropolitana de Lisboa:** a região de Lisboa em NUTS 2013. Em NUTS 2024 está dividida em `Grande Lisboa` e `Península de Setúbal`.
- **Oeste e Vale do Tejo:** região criada em NUTS 2024 com o Oeste, o Médio Tejo e a Lezíria do Tejo. Apesar do nome, não inclui Lisboa.
- **Agregação por municípios:** as regiões são sempre somadas a partir dos seus municípios, com a mesma lista aplicada a todos os anos. Mantém a série contínua, mas os totais não coincidem exactamente com os publicados pelo INE.
- **Revisão da população (2021):** o INE publica duas estimativas que não coincidem; a aplicação usa a revista a partir de 2021. Todas as taxas desde esse ano mudaram, e há um degrau de cerca de 1,7% em 2020/2021.
