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

**A lista de municípios de cada região é sempre a mesma em todos os anos.** Quando escolhe `Alentejo`, a aplicação usa os municípios que pertencem ao Alentejo na definição escolhida, e aplica essa mesma lista a todos os anos da série — é isso que a população usa sempre. Os óbitos vêm, por predefinição, das linhas regionais do INE quando existem (ver *Óbitos das Regiões*, abaixo). Isto tem uma razão e uma consequência.

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

### Óbitos das Regiões: Linhas do INE ou Soma dos Municípios

Ao lado da definição NUTS, no topo da página, há um segundo controlo:
`Óbitos das regiões`. Escolhe de onde vêm os óbitos de uma região.

**Porque existe.** O INE publica os óbitos por causa de cada município com os
totais completos, mas com a repartição por idade incompleta — sobretudo onde os
números são pequenos. Como todas as taxas da aplicação são calculadas por idade,
somar os municípios de uma região perde óbitos. No cancro do pulmão, em 2013, a
soma dos municípios fica 18% abaixo nos Açores, 12% na Madeira e 9% no Alentejo.
Em **2014** a perda chega a 30–84% em todas as regiões.

As linhas regionais publicadas pelo próprio INE não têm este problema.

| Opção | De onde vêm os óbitos | Quando usar |
|---|---|---|
| **Linhas regionais do INE** (predefinição) | a linha do INE para a região, sempre que existe; a soma dos municípios nos restantes anos | quase sempre |
| **Soma dos municípios** | sempre a soma dos municípios | para comparar com resultados antigos da aplicação, ou quando precisa de uma única fonte sem saltos |

A população é sempre a soma dos municípios, nas duas opções.

**Todas as regiões usam linhas do INE em todos os anos**, com uma excepção.
`Continente`, `Norte`, `Algarve`, `Açores` e `Madeira` são o mesmo território nas
duas definições e usam a sua própria linha. As regiões redesenhadas em 2024 —
`Centro`, `Alentejo`, `Oeste e Vale do Tejo` e a Área Metropolitana de Lisboa —
são **compostas a partir das sub-regiões (NUTS III)** nos anos em que a sua linha
não existe: por exemplo, o Alentejo de 2024 antes de 2022 é a soma das suas quatro
sub-regiões da definição de 2013. Cada composição foi verificada contra 2022, o
único ano publicado nas duas definições, e coincide exactamente.

A excepção são a `Grande Lisboa` e a `Península de Setúbal` antes de 2022: na
definição de 2013 a Área Metropolitana de Lisboa era uma única sub-região e não
pode ser dividida, por isso esses anos usam a soma dos municípios. Na prática a
diferença é mínima, porque os municípios de Lisboa são grandes, mas a aplicação
avisa.

**O que não fica corrigido.** Um município seleccionado sozinho não tem linha
mais fina onde ir buscar a idade. Quando escolhe municípios e uma causa
específica, a aplicação avisa que os valores podem estar subestimados, e que 2014
não é fiável ao nível municipal.

Na tabela de fontes, os óbitos que vieram de uma linha regional aparecem
identificados como `(linha regional)`, ou `(linhas regionais compostas)` quando a
região foi construída a partir de sub-regiões.

### ULS e ARS

Além das regiões NUTS, a lista de locais inclui a geografia do sistema de saúde:
as cinco **ARS** (Norte, Centro, Lisboa e Vale do Tejo, Alentejo, Algarve) e as
**ULS**. Tal como as regiões NUTS, são somadas a partir dos seus municípios e estão
disponíveis em todos os separadores — SMR, mortalidade evitável, mortalidade
infantil, previsões.

As ULS e as regiões NUTS **não coincidem**: cinco ULS atravessam uma fronteira NUTS
II, e a ARS Norte não é o mesmo território que a região NUTS Norte. Por isso aparecem
como entradas próprias, com o prefixo `ULS` ou `ARS`. A lista de municípios de cada
ULS é a mesma nas duas definições NUTS.

**Cinco ULS aparecem agrupadas.** Lisboa, Loures e Porto estão divididos entre duas
ULS ao nível da freguesia, e a aplicação não tem dados abaixo do município. Em vez de
as mostrar individualmente com valores errados, aparecem os agrupamentos exactos:

- `ULS Santo António + São João` — Gondomar, Maia, Porto, Valongo;
- `ULS Loures/Odivelas + Santa Maria + São José` — Lisboa, Loures, Mafra, Odivelas.

Estes valores não são comparáveis com os que o ficheiro de indicadores PNS2030
apresenta para cada uma destas ULS: esse ficheiro atribui o município partilhado
inteiro a cada ULS, e por isso conta Porto, Lisboa e Loures duas vezes.

**Óbitos por causa nas ULS.** Sete ULS coincidem exactamente com uma sub-região NUTS
III (Alto Minho, Viseu Dão-Lafões, Alentejo Litoral, Baixo Alentejo, Alto Alentejo,
Alentejo Central e Algarve) e usam as linhas regionais do INE. As restantes são
somas de municípios e têm o mesmo cuidado dos municípios isolados — a aplicação
avisa. Em 2014 esse cuidado é sério: a ULS Guarda, por exemplo, aparece com zero
óbitos por cancro do pulmão nesse ano.

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
- agregar vários anos ajuda, mas não cria acontecimentos que não houve;
- **em 1995-2001 as regiões, ULS e municípios estão subestimados.** Antes de 2011 o
  INE não publica contagens completas de óbitos com menos de 1 ano por município;
  a aplicação usa então a idade «menos de 1 ano» das estatísticas de causas de
  morte, cujos municípios somam cerca de 85% do total nacional nesses anos (e
  praticamente 100% em 2002-2010). A aplicação avisa quando a selecção inclui
  esses anos, e o separador de Indicadores de Planeamento marca o valor com `†`.
  Portugal lê a sua própria linha e não é afectado.

**Correcção de Setembro de 2026.** Até esta versão, os nados-vivos de 1995-2013
estavam errados para Lisboa: a linha «Lisboa» do INE nesses anos é a região, não o
município, e a aplicação registava cerca de 37.000 nascimentos em vez de cerca de
5.600. A taxa de Lisboa, e de qualquer região ou ULS que a incluísse, ficava muito
abaixo da real. Os nados-vivos e os óbitos com menos de 1 ano são agora
identificados pelo código do município e não pelo nome, o que também separa as duas
Calhetas e as duas Lagoas. Desde 2011, os óbitos com menos de 1 ano vêm das
contagens completas do INE, e já não da repartição por idade das causas de morte
(que em 2014 perdia metade dos óbitos municipais).

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

## 10-A. Separador Mortalidade Evitável

Mortalidade evitável são os óbitos **antes dos 75 anos** que poderiam ter sido
evitados. A distinção que este separador faz é a que o Eurostat e a OCDE usam
desde a revisão conjunta de 2019, e separa duas responsabilidades diferentes:

- **Prevenível** — evitável por saúde pública e prevenção primária: tabaco,
  álcool, segurança rodoviária, prevenção do suicídio. É matéria de política de
  saúde pública.
- **Tratável** — evitável por cuidados de saúde atempados e eficazes: rastreio,
  diagnóstico, tratamento. É matéria de organização e qualidade dos serviços.

Ler as duas em conjunto responde a perguntas diferentes, e é por isso que
aparecem separadas em vez de somadas num único indicador de "evitável".

Os 75 anos fazem parte da definição, não são uma opção: o separador aplica
sempre esse limite.

### O que a tabela mostra

Uma repartição dos óbitos com menos de 75 anos que **fecha no total**:

| Linha | O que é |
|---|---|
| Prevenível | as 16 causas atribuídas à prevenção primária |
| Tratável | as 17 causas atribuídas aos cuidados de saúde |
| Evitável (total) | a soma das duas |
| Não classificada (por resolver) `*` | seis causas deixadas de fora, ver abaixo |
| Não evitável / fora das listas | tudo o resto |
| Todas as causas, < 75 anos | o total |

Para cada linha são apresentados os óbitos, a percentagem do total com menos de
75 anos, e a taxa padronizada por 100.000 nessa faixa etária. A segunda tabela
lista exactamente que causas entram em cada grupo, e pode ser exportada.

### As seis causas excluídas

O Eurostat define as listas ao nível da CID-10. Os dados do INE vêm na lista
sucinta europeia, que é menos detalhada, e há seis rubricas que não podem ser
atribuídas a um dos grupos sem um juízo clínico:

| Causa | Porquê |
|---|---|
| Doenças isquémicas do coração | o Eurostat reparte-a em 50% prevenível e 50% tratável; a lista sucinta não permite separar |
| Doenças cérebro-vasculares | classificada como tratável, mas parte é atribuível a prevenção primária |
| Tumor maligno do estômago | consta de umas revisões da lista e não de outras |
| Tumor maligno do tecido linfático e hematopoético | só o linfoma de Hodgkin e a leucemia infantil são tratáveis; a rubrica agrega tudo |
| Tumor maligno do ovário | consta de umas revisões e não de outras |
| Hepatite viral | a hepatite B é prevenível por vacinação, a C é tratável |

**Estas seis causas ficaram de fora dos totais**, e não são uma nota de rodapé:
representam cerca de 18% dos óbitos com menos de 75 anos, dos quais dois terços
são as duas rubricas cardiovasculares. Por isso aparecem como linha própria na
tabela, marcada com `*`, em vez de desaparecerem.

### Como ler estes números

**São um limite inferior, não uma estimativa.** Além das seis excluídas, a lista
sucinta não tem rubrica para alguns grupos CID que constam das listas do
Eurostat — os 16 tumores que nomeia somam menos do que o seu próprio total de
tumores malignos. O que falta não está contabilizado.

**Servem para comparar, não para citar.** A base é consistente entre locais e
entre anos, pelo que as comparações são válidas. Mas os valores não reproduzem
os números publicados pelo Eurostat para Portugal e não devem ser apresentados
como tal.

Cuidados:

- em concelhos pequenos aplicam-se os mesmos cuidados de sempre: use a agregação
  plurianual e leia os intervalos;
- a definição NUTS activa, no topo da página, muda o que uma região significa;
- as taxas exigem denominador, pelo que os anos sem população publicada não
  estão disponíveis neste separador.

## 10-B. Separador Indicadores de Planeamento

Os indicadores demográficos e de mortalidade que os Planos Locais de Saúde
apresentam para cada ULS, calculados para qualquer local da aplicação: Portugal,
Continente, regiões NUTS, ARS, ULS e municípios. Seguem os indicadores do
ficheiro de apoio aos PLS (a referência `[I…]` de cada linha), recalculados a
partir dos dados actuais do INE.

O separador é construído à volta de **um local**. Escolhe-se o local, o indicador
e o intervalo de anos; o separador mostra esse indicador em gráfico e tabela, ao
lado das áreas que contêm o local, e resume todos os outros.

| Indicador | Ref. | Cálculo |
|---|---|---|
| População residente | I1 | estimativa anual do INE |
| Proporção de jovens, de idosos e de 75+ | I1 | grupo etário sobre a população total, em % |
| Índice de envelhecimento | I4 | 65 e mais anos por 100 com 0-14 anos |
| Índice de dependência de jovens | I5 | 0-14 anos por 100 com 15-64 anos |
| Índice de dependência de idosos | I6 | 65 e mais anos por 100 com 15-64 anos |
| População residente nos Censos e variação desde o censo anterior | I2 | Censos de 1991, 2001, 2011 e 2021 |
| População por nível de escolaridade mais elevado completo | I24 | % da população nos Censos: sem nível completo, básico, secundário, superior |
| Taxa de analfabetismo | I26 | % da população com 10 e mais anos que não sabe ler nem escrever, nos Censos |
| Nados-vivos | I7 | contagem |
| Taxa bruta de natalidade | I8 | nados-vivos por 1.000 habitantes |
| Óbitos | I37 | contagem, todas as causas e idades |
| Taxa bruta de mortalidade | I38 | óbitos por 1.000 habitantes |
| Índice sintético de fecundidade | I9 | soma das taxas de fecundidade por idade da mãe (15-49), vezes 5 |
| Nascimentos em mães com menos de 20 anos | I32 | % dos nados-vivos, no triénio |
| Nascimentos em mães com 35 e mais anos | I33 | % dos nados-vivos, no triénio |
| Nascimentos pré-termo | I35 | menos de 37 semanas, % dos nascimentos com duração conhecida, no triénio |
| Nascimentos com baixo peso | I36 | menos de 2.500 g, % dos nascimentos com peso conhecido, no triénio |
| Beneficiários do RSI | I13, I14 | contagem, e por 1.000 habitantes com 15 e mais anos |
| Pensionistas da segurança social | I15, I16 | contagem, e por 1.000 habitantes com 15 e mais anos |
| Valor médio das pensões | I17 | valor total das pensões sobre o total de pensionistas (€ por ano) |
| Poder de compra per capita | I28 | Portugal = 100; anos ímpares (estudo bienal) |
| Ganho médio mensal | I27 | ganho médio dos trabalhadores por conta de outrem, ponderado pelo número de trabalhadores (€/mês) |
| Trabalhadores por conta de outrem e repartição por sector | I12 | contagem e % nos sectores primário, secundário e terciário |
| Resíduos urbanos por habitante | I64, I65 | total e recolha selectiva, kg por habitante |
| Esperança de vida à nascença e aos 65 anos, total e por sexo | I10 | tábua de mortalidade abreviada por triénio (ver abaixo) |
| Taxa de mortalidade infantil | I39 | óbitos com menos de 1 ano por 1.000 nados-vivos, no triénio |
| Mortalidade neonatal, neonatal precoce e pós-neonatal | I40-I42 | óbitos com menos de 28 dias, menos de 7 dias e 28-364 dias por 1.000 nados-vivos, no triénio |
| Mortalidade fetal tardia | I43 | fetos-mortos com 28 ou mais semanas por 1.000 nascimentos (nados-vivos + fetos-mortos), no triénio |
| Mortalidade perinatal | I44 | fetos-mortos com 28 ou mais semanas e óbitos com menos de 7 dias, por 1.000 nascimentos, no triénio |
| Mortalidade proporcional por grandes grupos de causas | I45, I46 | óbitos de cada grupo sobre o total, no triénio, para todas as idades e para as idades abaixo de 75 |

Os indicadores aparecem agrupados por tema: Demografia, Natalidade, Contexto
social, Ambiente e Mortalidade.

Ficam de fora, por agora: os indicadores que não vêm do INE (I11 do IEFP, I18 e
I19 da PORDATA, I34 e I49 do SIM@SNS), o abastecimento de água e a drenagem de
águas residuais (I61 e I62), cuja série municipal do INE terminou em 2009, e a
água segura (I63), publicada como percentagem sem denominador, que por isso não
se pode agregar. As folhas I29, I30 e I31 estão marcadas como descontinuadas no
próprio ficheiro.

Os indicadores dos Censos (I2, I24, I26) existem apenas nos anos censitários
(1991, 2001, 2011 e 2021).

### Comparadores

Para cada local, a aplicação propõe as áreas que o contêm, uma por nível:

| Local escolhido | Comparadores propostos |
|---|---|
| Município | a sua ULS, ARS, NUTS III, NUTS II, NUTS I e Portugal |
| ULS | a ARS, a NUTS III ou NUTS II em que cabe inteira, o Continente e Portugal |
| ARS ou NUTS II | o NUTS I (Continente ou região autónoma) e Portugal |
| NUTS I | Portugal |

Todos começam seleccionados; pode desmarcar os que não interessam. Uma área com
exactamente os mesmos municípios do local não é proposta, porque repetiria os
mesmos valores: o município de Matosinhos não tem a ULS Matosinhos como
comparador. As regiões autónomas não têm ULS nem ARS na aplicação.

**Só os indicadores comparáveis têm comparadores**: taxas, proporções, índices e
valores por habitante, que não dependem do tamanho da área. As contagens
(população, nados-vivos, óbitos, beneficiários do RSI, pensionistas) mostram-se
apenas para o local, porque o número de um município não diz nada ao lado do da
sua região. Ao escolher uma contagem, os comparadores desaparecem e o gráfico
passa a barras.

Cada nível tem sempre a mesma cor em todos os gráficos: o local a azul, a ULS a
laranja, a ARS a verde-água, a NUTS III a amarelo, a NUTS II a rosa, o NUTS I a
verde e Portugal a violeta.

### Portugal: total do INE ou soma dos municípios

A opção **Portugal** escolhe o que representa Portugal em todo o separador
(comparadores, significância, classificação das ULS, funil e perfil):

- **Total publicado pelo INE**: inclui os acontecimentos de residentes cujo
  município é desconhecido (0,3% a 0,9% dos óbitos).
- **Soma dos 308 municípios**: sem esses acontecimentos, como qualquer região,
  ULS ou município. Compara igual com igual, e é a escolha certa quando se
  pergunta se um local se afasta do país.

A diferença é pequena (óbitos de 2023: 118.344 no total do INE, 118.328 na soma
dos municípios; a taxa bruta muda na terceira casa decimal), mas pode decidir um
teste de significância no limite, sobretudo para as regiões grandes. O Excel de todas as áreas traz as duas linhas.

### Significância face a Portugal

Nas tabelas, cada valor com intervalo de confiança leva uma marca:

| Marca | Significado |
|---|---|
| ▲ | o intervalo de 95% fica inteiramente **acima** do valor de Portugal |
| ▼ | o intervalo de 95% fica inteiramente **abaixo** |
| = | o intervalo inclui o valor de Portugal: sem diferença demonstrável |

É o critério do PHE Fingertips. Só se aplica a taxas, proporções e esperança de
vida; as contagens e os índices sem intervalo não levam marca. A marca não diz
se a diferença é boa ou má: uma esperança de vida acima é boa, uma mortalidade
acima não. Na classificação das ULS as barras têm a mesma leitura: laranja acima,
azul abaixo, cinzento sem diferença; a ULS do local tem contorno escuro.

### Funil

O subseparador **Funil** mostra todas as ULS, municípios ou NUTS III de uma vez:
cada ponto é o valor contra o tamanho do denominador (nados-vivos, população). As
linhas marcam o que uma unidade desse tamanho mostraria **só por acaso** à volta
de Portugal: 95% (tracejado) e 99,8% (pontilhado). As unidades pequenas espalham-se
muito dentro das linhas, as grandes pouco; um ponto fora das linhas afasta-se
mais do que o acaso explica. Com 308 municípios, cerca de 15 ficariam fora do
limite de 95% por acaso; fora do de 99,8%, menos de um. Um losango marca os que
estão fora do limite de 99,8%.

Os limites são os quantis exactos da contagem sob a taxa de Portugal (Poisson
para taxas de acontecimentos, binomial para proporções), interpolados
(Spiegelhalter, 2005). Assim, um município pequeno sem nenhum óbito infantil
nunca aparece como «significativamente baixo». O funil só existe para
indicadores com um modelo de contagem; para as taxas brutas, a dispersão reflecte
sobretudo a estrutura etária.

### Escolaridade por idade

A opção **Escolaridade [I24], população** restringe os quatro indicadores de
escolaridade a quem tem uma idade mínima (15, 20, 25… 75 e mais anos). Com toda
a população, como no ficheiro de apoio, as crianças contam como «sem nível de
escolaridade completo», o que faz parecer menos escolarizadas as áreas com mais
crianças. Com 25 e mais anos, por exemplo, compara-se a escolaridade de adultos.

O INE só publica a escolaridade por idade e município nos Censos de 2011 e 2021;
com idade mínima, 1991 e 2001 ficam sem valor. A escolha vale para o gráfico, o
resumo, o Excel do local e o perfil.

### Notas com asterisco

Um `*` no fim do nome de um indicador remete para uma nota, mostrada por baixo do
gráfico, no fim do resumo, no Excel e no perfil:

- **Ganho médio e trabalhadores por sector [I27, I12]**: Quadros de Pessoal,
  contados no local de trabalho e não no de residência, sem a Administração
  Pública nem os trabalhadores por conta própria.
- **Esperança de vida [I10]**: reproduz o Eurostat e fica cerca de 0,8-0,9 anos
  acima dos valores do INE, que usa outra metodologia; compare valores da
  aplicação entre si.

### Mortalidade padronizada, prematura e evitável

O tema **Mortalidade padronizada** acrescenta indicadores que não estão no
ficheiro de apoio, todos por triénio:

| Indicador | O que mede |
|---|---|
| Razão padronizada de mortalidade (SMR) | óbitos observados sobre os esperados com as taxas por idade de Portugal no mesmo triénio; Portugal = 100 |
| Taxa de mortalidade padronizada | a taxa que o local teria com a População Padrão Europeia de 2013, por 100.000 habitantes |
| Mortalidade prematura padronizada | o mesmo, só antes dos 75 anos |
| Mortalidade evitável por prevenção / por cuidados de saúde | taxas padronizadas antes dos 75 anos, com as listas Eurostat/OCDE de 2019 |
| Óbitos prematuros e óbitos evitáveis | contagens antes dos 75 anos |
| Anos potenciais de vida perdidos | anos que faltavam até aos 70 em cada óbito, por 100.000 residentes com menos de 70 anos |

Ao contrário da taxa bruta, estes indicadores comparam áreas com estruturas
etárias diferentes: um concelho envelhecido tem mais óbitos sem ter por isso
mais mortalidade. A SMR segue a opção **Portugal** (total do INE ou soma dos
municípios).

O subseparador **Mortalidade por causa** mostra, para os 13 grandes grupos de
causas e para todas as causas, os óbitos observados e esperados, a SMR com o
seu intervalo e as taxas padronizadas (todas as idades e antes dos 75), para
ambos os sexos, homens ou mulheres.

**Como se lidam os óbitos sem idade.** O INE publica os óbitos por idade e causa
de cada município incompletos: em 2023 os grupos etários municipais têm 97,9%
dos óbitos por doenças do aparelho circulatório, e em 2014 só 65% dos
suicídios; o total de todas as idades, esse, é completo. Somados tal como
estão, todas as regiões e ULS ficariam 1% a 3% abaixo de Portugal. Por isso os
óbitos de cada município e causa são primeiro completados até ao total, e os
que faltam distribuídos pelas idades com o perfil dos óbitos que faltam no
país (a linha de Portugal por idade menos a soma dos municípios). Assim a soma
dos municípios reproduz Portugal por idade, e o Norte reproduz a linha do INE
de óbitos antes dos 75 anos (33.155 em 2022-2024). A marca `‡` assinala os
valores em que mais de 2% dos óbitos foram redistribuídos.

**Mortalidade evitável é um limite inferior**: seis causas da lista sucinta
(cerca de 18% dos óbitos antes dos 75) não se podem atribuir a um grupo sem
juízo clínico e ficam de fora. Serve para comparar locais e anos, não para
reproduzir os números do Eurostat.

O funil também existe para a SMR, com os óbitos esperados como dimensão.

### Cuidados de saúde primários (Portal da Transparência do SNS)

O subseparador **Cuidados de saúde primários** mostra, por ULS, indicadores
mensais do Portal da Transparência do SNS: utentes sem médico de família,
utilização de consultas, rastreios oncológicos (mama, colo do útero, cólon e
reto), diabetes (exame dos pés, HbA1c), hipertensão e vigilância do
recém-nascido. Para o local escolhido mostra a sua ULS (um município mostra a
ULS a que pertence), a ARS e o Continente, e em baixo todas as ULS no último
período completo, com a significância face ao Continente.

Três cuidados na leitura:

- **A base são os utentes inscritos** nos cuidados de saúde primários, não os
  residentes.
- **Vários indicadores acumulam ao longo do ano ou do semestre** e recomeçam: os
  rastreios e o exame dos pés sobem de Janeiro a Dezembro; a tensão arterial e a
  HbA1c de Janeiro a Junho e de Julho a Dezembro. Um valor de Julho não se
  compara com um de Dezembro. O gráfico parte a linha em cada recomeço e marca
  os fins de ciclo, os únicos valores comparáveis; a tabela compara o ciclo em
  curso com o mesmo mês do ano anterior.
- **Só há ULS desde Janeiro de 2024** (antes, o portal publica por ACES), nada
  para as regiões autónomas, e o último mês é provisório. As cinco ULS de Lisboa
  e do Porto somam-se nos dois agrupamentos exactos da aplicação. Uma NUTS III
  que não coincide com um conjunto de ULS não tem valores.

### Mortalidade semanal e excesso de mortalidade

O subseparador **Mortalidade semanal** mostra os óbitos semanais publicados pelo
INE, com poucas semanas de atraso, para a NUTS III do local (o INE não os
publica por município nem por ULS; uma área que atravessa várias NUTS III
mostra a região maior que a contém). A linha azul são os óbitos observados; a
faixa, os esperados com o seu intervalo de 95%.

Os **esperados** aplicam à população do ano as taxas de mortalidade semanais por
idade (menos de 65, 65-74, 75-84, 85 e mais) dos anos de base. Assim o
envelhecimento e o crescimento da população não aparecem como excesso. Os anos
de base são até cinco anteriores, **desde 2023**: depois do excesso da COVID-19
(2020-2022) e na mesma série de população (o INE reviu a população a partir de
2021; com as taxas de 2018-2019 sobre a população antiga, todos os anos
recentes pareciam ter menos óbitos do que o esperado). Por isso 2024 ainda não
tem esperados, 2025 tem dois anos de base e 2026 três.

O gráfico de baixo acumula o excesso ao longo de cada ano; a tabela dá o total
do ano até à última semana, com o intervalo. As últimas semanas são provisórias
e sobem à medida que chegam os registos em atraso.

### Actualização automática

Os dados que mudam com frequência (óbitos semanais, Portal do SNS, e os óbitos,
a população e os óbitos com menos de 1 ano mais recentes do INE) são
actualizados todas as segundas-feiras por uma tarefa agendada no computador
onde a aplicação corre (o INE recusa ligações dos servidores do GitHub). Cada
actualização fica registada no histórico dos dados com a sua data, e as versões
substituídas são guardadas.

### Perfil do local (Word)

O botão **Perfil do local (Word)** gera um documento editável com: um resumo dos
indicadores acima e abaixo de Portugal, o quadro de todos os indicadores com os
comparadores e as marcas de significância, a pirâmide etária, a evolução de seis
indicadores, a posição da ULS do local entre as ULS do Continente, a mortalidade
proporcional e as notas de método. Usa as mesmas opções do separador (local,
comparadores seleccionados, Portugal, escolaridade, último ano). Demora cerca de
20 segundos.

### Os subseparadores

- **Indicador** — o gráfico do indicador ao longo dos anos escolhidos e, por
  baixo, a tabela com os valores e os intervalos de confiança (o período mais
  recente primeiro). Para indicadores comparáveis, linhas: o local com a sua faixa
  de intervalo de confiança, e uma linha por comparador. Para contagens, barras do
  local com o intervalo como barra de erro. Pontos vazios têm uma marca (`*` ou
  `†`); uma linha pontilhada vertical marca uma mudança de série.
- **Resumo do local** — todos os indicadores no último período disponível até ao
  último ano escolhido, com os comparadores ao lado dos que são comparáveis.
- **Comparação entre ULS** — as 34 ULS do Continente e os dois agrupamentos
  exactos, ordenados, com Portugal como linha tracejada. A ULS do local (ou o
  próprio local, se for uma ULS) aparece a azul. Só para indicadores comparáveis.
- **Pirâmide etária** — a estrutura por idade e sexo do local, em percentagem da
  sua população, com o contorno do primeiro comparador seleccionado (ou de
  Portugal) por cima, para comparar estruturas de áreas de tamanhos diferentes (I3).
- **Mortalidade proporcional** — os 13 grandes grupos de causas, mais «Restantes
  causas» para que o total feche em 100%: barras para o local, pontos para os
  comparadores, e a tabela por baixo. Pode escolher todas as idades [I45] ou
  apenas os óbitos antes dos 75 anos [I46].
- **Notas** — definições, fontes, mudanças de série e diferenças face ao ficheiro
  de apoio.

Todos os gráficos mostram o valor exacto ao passar o rato sobre um ponto ou barra.

### Descarregar em Excel

No fim do separador há dois ficheiros Excel:

| Ficheiro | Conteúdo |
|---|---|
| **Local e comparadores** | o local e os comparadores seleccionados, nos anos escolhidos: folha Leia-me, Resumo, uma folha por indicador (áreas em linhas, anos em colunas), Dados em formato longo com intervalos de confiança, numerador e denominador, Pirâmide etária e Mortalidade proporcional |
| **Todas as áreas** | Portugal, NUTS I, II e III, as 5 ARS, as 36 ULS e agrupamentos e os 308 municípios, todos os indicadores e anos, uma folha por indicador, mais a pirâmide (2011, 2021 e o último ano) e a mortalidade proporcional de todos os triénios. É o equivalente automático do ficheiro de apoio aos PLS |

Nos dois, a folha Leia-me indica a **data de importação dos dados** e a definição
das regiões usada, e os valores com uma marca aparecem a cinzento e itálico. O
ficheiro de todas as áreas é gerado da primeira vez que é pedido (cerca de meio
minuto) e guardado até à importação seguinte.

### Como são construídos os valores

Cada local é a soma dos seus municípios, e cada indicador é uma razão dessas
somas, nunca uma média de taxas municipais. Portugal e o Continente usam as
linhas publicadas pelo INE, que incluem os acontecimentos cujo município de
residência é desconhecido. Por isso a soma das regiões fica ligeiramente abaixo
de Portugal, cerca de 0,3% a 0,9% dos óbitos.

Os óbitos deste separador vêm do **total de todas as idades** que o INE publica
por município, e não da repartição por idade usada nos outros separadores. Essa
repartição está incompleta ao nível municipal, sobretudo em 2014. Contagens,
taxas brutas e proporções de todas as idades não precisam da idade, pelo que
aqui são exactas em qualquer nível.

### População média e lacunas nos dados do INE

Nascimentos, óbitos, beneficiários do RSI e resíduos acontecem ao longo do ano,
por isso as suas taxas dividem pela **população média** do ano (a média das
estimativas a 31 de Dezembro do ano anterior e do próprio ano), como o INE e o
ficheiro de apoio. Com isso as taxas brutas de natalidade e mortalidade e os
resíduos por habitante de 2024 coincidem com os do INE em todos os municípios.
Os pensionistas, contados a 31 de Dezembro, dividem pela população nessa data.

Uma célula em branco no INE não é um zero:

- Nos resíduos, pensões, trabalhadores e ganho médio, publicados para todos os
  municípios, um município sem valor deixa **sem valor** as áreas que o contêm
  (por exemplo, os Açores em 2013-2014 nos trabalhadores por sector).
- **Odivelas, Trofa e Vizela** foram criados em 1998. Até esse ano os seus
  nascimentos e óbitos estão registados em Loures, Santo Tirso e Guimarães; só
  há valores para áreas com os dois municípios de cada par, também na esperança
  de vida e na mortalidade proporcional dos triénios que incluem esses anos. Os
  resíduos de Loures incluem sempre os de Odivelas (serviço conjunto SIMAR).
- Quando o INE deixa em branco o total de óbitos de um município (Vimioso em
  2024), usa-se a soma dos grupos etários publicados.
- Nos **trabalhadores por sector**, o INE oculta dois sectores quando um deles
  revelaria uma empresa. O remanescente é repartido na proporção do resto da
  NUTS III, e a marca `≈` assinala quotas com mais de 1% de trabalhadores
  estimados.

### Intervalos e asterisco

Contagens e taxas de acontecimentos (nascimentos, óbitos, mortalidade infantil,
proporções) têm intervalo de confiança de 95%. Os índices de estrutura da
população não têm: as estimativas de população não são uma amostra de
acontecimentos. O `*` na mortalidade infantil marca um triénio com menos de
1.000 nados-vivos, como no resto da aplicação.

### Ganho médio mensal e trabalhadores por sector

Vêm dos Quadros de Pessoal (MTSSS/GEP), publicados pelo INE por município, e
cobrem os trabalhadores **por conta de outrem**, contados no **local de
trabalho** e não no local de residência. Os indicadores do ficheiro de apoio
usam, para o I12, a população empregada dos Censos, por residência: as duas
medidas não são comparáveis entre si, embora respondam à mesma pergunta.

O ganho médio de uma área é a média dos valores municipais **ponderada pelo
número de trabalhadores**, não a média simples: um município pequeno não pesa o
mesmo que a sede de concelho onde estão os empregos. Para uma ULS com um só
município o valor é o do próprio município, e coincide exactamente com o do
ficheiro de apoio (ULS Matosinhos, 2013-2018).

Cobertura: ganho médio desde 2011, trabalhadores por sector desde 2013. O ganho
médio só é apresentado a partir de 2013, porque antes disso não há os pesos.

### Indicadores dos Censos

Vêm das séries históricas do INE por município (população residente desde 1864,
escolaridade desde 1940, alfabetismo desde 1878), pelo que os valores dos Censos
de 1991, 2001, 2011 e 2021 são comparáveis entre si.

A taxa de analfabetismo segue a definição do INE — população com 10 e mais anos
que não sabe ler nem escrever, sobre a população com 10 e mais anos. O INE publica
a taxa por município, não a contagem; para um agrupamento, a aplicação faz a média
das taxas dos municípios ponderada pela população com 10 e mais anos de cada um,
o que é equivalente à razão das contagens.

Na escolaridade, a série conta apenas quem tem um nível completo; «sem nível
completo» é a diferença para a população total do censo.

### Mortalidade proporcional abaixo dos 75 anos

A repartição por todas as idades [I45] lê os totais de óbitos por município e
causa, que estão completos. A repartição abaixo dos 75 anos [I46] precisa da
idade dos óbitos, e aí cada área usa a melhor fonte disponível:

- **a linha regional do INE**, onde existe: Portugal, Continente, todas as regiões
  NUTS e as ULS que coincidem com uma NUTS III (Alto Minho, Viseu Dão-Lafões,
  Algarve e as quatro do Alentejo);
- **a soma dos municípios**, nas restantes áreas.

A soma dos municípios é fiável: em 2020-2022 reproduz exactamente os óbitos abaixo
dos 75 anos do Alto Minho e do Algarve publicados pelo INE, e as quotas com um
desvio máximo de 0,3 pontos percentuais.

A excepção é 2014, ano em que o INE publicou a idade de apenas 80% dos óbitos por
município (52% no pior grupo de causas). Nos três triénios que incluem 2014
(2012-2014, 2013-2015 e 2014-2016), as áreas sem linha regional aparecem marcadas
com `§`: as suas quotas podem estar desviadas em cerca de 2 pontos percentuais.

Os óbitos não são reescalados para os totais completos. Testado contra as linhas
do INE, reescalar piora as quotas em vez de as melhorar: os óbitos sem idade
publicada estão sobretudo nas idades mais altas, pelo que distribuí-los na
proporção dos restantes coloca demasiados abaixo dos 75 anos.

### Mortalidade fetal tardia e perinatal

O INE publica por município os óbitos perinatais (fetos-mortos com 28 ou mais
semanas de gestação mais óbitos com menos de 7 dias), mas não os fetos-mortos em
separado. A aplicação obtém-nos subtraindo aos óbitos perinatais os óbitos com
menos de 7 dias do mesmo município, que conhece desde 2011. Por isso estes dois
indicadores começam no triénio 2011-2013.

### Esperança de vida à nascença

Apresentada à nascença e aos 65 anos, no total e por sexo, porque os planos usam
as duas: a primeira resume a mortalidade em todas as idades, a segunda a
mortalidade depois dos 65.

Calculada com uma tábua de mortalidade abreviada (método de Chiang), por triénio,
a partir dos óbitos de todas as causas por grupo quinquenal de idade até «85 e mais
anos» e da população a meio do ano. É o método usado pelo Eurostat e pela Public
Health England, e tem intervalo de confiança de 95%. Não é apresentada para áreas
com 5.000 habitantes-ano ou menos no triénio, nem quando o intervalo excede 20
anos.

Alguns municípios não têm todos os óbitos repartidos por idade no INE (sobretudo
em 2014). Esses óbitos são distribuídos pelas idades na proporção dos restantes do
mesmo município, e o valor é marcado com `‡` quando representam mais de 2% dos
óbitos do triénio.

**Não compare estes valores com os publicados pelo INE.** Para Portugal a aplicação
coincide com o Eurostat, à nascença e aos 65 anos (2017-2019: 81,9 e 20,7 anos;
Eurostat 82,0 e 20,6 em 2019), mas o INE,
com a sua Metodologia 2007, publica valores cerca de 0,8-0,9 anos mais baixos. A
diferença é praticamente constante entre regiões: nas 26 NUTS III de 2021-2023 a
correlação com o INE é 0,97. Os valores da aplicação servem para comparar áreas e
anos entre si.

### Porque podem diferir do ficheiro de apoio aos PLS

- **População revista.** Desde 2021 a aplicação usa a série revista do INE.
  Valores calculados com a estimativa anterior ficam desactualizados: o índice
  de envelhecimento do Alto Minho em 2024 era 270,3 com a estimativa antiga e é
  240,2 com a revista.
- **ULS que partilham um município.** Lisboa, Loures e Porto estão divididos
  entre ULS ao nível da freguesia. Atribuir o município inteiro a cada ULS conta
  a mesma população duas vezes: no ficheiro, a soma das ULS do Norte
  ultrapassa a ARS Norte em cerca de 3.000 óbitos. A aplicação mostra antes os
  dois agrupamentos exactos.
- **Ganho médio mensal das ULS com vários municípios.** No ficheiro, a linha de
  uma ULS repete o valor do primeiro município por ordem alfabética: a série da
  ULS Alto Minho de 2013 a 2018 (802,3; 797,3; 801,3; 815,8; 875,8; 882,1) é,
  ano a ano, a de Arcos de Valdevez. A aplicação pondera pelos trabalhadores
  (878,9 em 2013). Onde a ULS tem um só município, os valores coincidem.
- **Poder de compra de uma ULS.** É a soma das quotas dos municípios no poder de
  compra nacional sobre a soma das suas quotas de população, e não a média dos
  índices municipais. Nalgumas ULS o ficheiro tem valores que não coincidem com os
  do INE: Matosinhos em 2021 tem 118,1 no INE e 130,6 no ficheiro.
- **Mortalidade neonatal das ULS.** Na aplicação, a neonatal e a pós-neonatal somam
  sempre a infantil. No ficheiro isso não acontece em várias ULS (Alto Minho,
  2022-2024: 1,1 + 0,9 contra 2,4), o que sugere linhas desalinhadas.
- **Índice sintético de fecundidade.** O do ficheiro para 2024 (Continente 1,41)
  é anterior à revisão da população; o INE publica agora 1,27 para Portugal, valor
  que a aplicação reproduz.
- **Mães com menos de 20 anos.** O ficheiro conta só as mães de 15-19 anos; a
  aplicação inclui também as de 10-14, como diz o nome do indicador (Continente
  2022-2024: 1,87% contra 1,8%).
- **Pensões.** Os valores coincidem com o ficheiro (Série 2017 da segurança
  social). Em 2017 há uma mudança de série, com cerca de 5,5% menos
  pensionistas; o gráfico de evolução marca-a.
- **Edição dos dados de 2022.** A aplicação lê 2022 da tabela NUTS 2024 do
  INE. Para o Continente em 2020-2022 obtém 356.355 óbitos, contra 356.333 no
  ficheiro: uma diferença de 22 (0,006%), e nula em vários grupos de causas
  (circulatório, geniturinário, perinatal).

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

### Histórico dos dados

O INE revê dados que já publicou: os anos mais recentes começam por ser
provisórios, as estimativas de população são re-estimadas (como em 2021), e
indicadores são substituídos por novas edições. Por isso a mesma análise, feita
hoje e daqui a um ano, pode dar valores diferentes, sem que nada na análise tenha
mudado.

A aplicação guarda o que é preciso para o explicar:

- **A data de importação dos dados** aparece no topo de todas as páginas («Dados
  importados até …»). Ao guardar ou partilhar um resultado, anote essa data.
- **Cada importação fica registada**: que ficheiros acrescentou, quais reviu, e
  em quanto mudaram (linhas com valor alterado e o total de Portugal antes e
  depois).
- **A versão anterior de qualquer ficheiro revisto é guardada**, e não
  substituída. Voltar a importar um ano que o INE não reviu não altera nada.

No fim deste separador, a secção «Histórico dos dados» mostra:

| Quadro | O que diz |
|---|---|
| Conjuntos de dados | para cada conjunto, os anos, a primeira e a última importação, e a última vez que valores já existentes mudaram |
| Importações | cada importação, com os ficheiros novos e revistos |
| Valores revistos | para a importação escolhida, por conjunto e ano: linhas alteradas e o total de Portugal antes e depois |

O total de Portugal pode não mudar quando a correcção é entre municípios. Foi o
caso da correcção dos nados-vivos de Lisboa em Setembro de 2026: Portugal manteve
o total, mas as linhas alteradas mostram a revisão.

**Repetir uma análise com os dados de uma data anterior.** Quem tenha a cópia do
repositório pode reconstruir os dados tal como estavam num dia e abrir a
aplicação sobre eles:

```sh
Rscript tools/data_as_of.R date=2026-09-01
MORTALITY_SNAPSHOT_DIR=.mortality-shiny-cache/data_as_of/2026-09-01/snapshots Rscript -e 'shiny::runApp()'
```

A reconstrução ocupa pouco espaço: os ficheiros que não mudaram não são
copiados. Os mapas de municípios por região e por ULS são sempre os actuais.

## 12. Exportação de Resultados

Todas as tabelas podem ser descarregadas em CSV e os gráficos estáticos em PNG. O
separador de Indicadores de Planeamento exporta ainda dois ficheiros Excel (ver a
secção 10-B).

**Todos os ficheiros dizem de que versão dos dados vieram.** No fim de cada CSV há
uma linha de comentário, e no rodapé de cada PNG uma legenda, com a data de
importação dos dados do INE, a definição das regiões em uso e a data de
exportação:

```text
# Dados do INE importados até 2026-09-17; regiões NUTS 2024; exportado em 2026-09-18
```

A linha começa por `#` e fica depois dos dados, pelo que o ficheiro continua a
abrir normalmente no Excel; em R ou Python pode ser ignorada com a opção de
comentário (`read.csv(..., comment.char = "#")`). Guarde-a junto dos resultados:
é o que permite explicar, mais tarde, diferenças entre análises feitas em datas
diferentes (ver a secção sobre o histórico dos dados, no separador Disponibilidade
de Dados).

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

### As Mudanças de Base da População

A população do arquivo vem de três indicadores do INE, encadeados. Onde um
termina e o outro começa há uma **mudança de base**, e os dois não coincidem
exactamente nos anos que ambos publicam. Há duas, em 2013/2014 e em 2020/2021.

**2013/2014.** Nos anos que ambos publicam, os totais de população coincidem a
menos de 0,2%. Mas as **taxas padronizadas** diferem cerca de 3%, porque os dois
indicadores distribuem a população pelas idades de forma diferente e a
padronização é sensível a isso. O efeito é visível:

| Portugal, taxa padronizada | |
|---|---:|
| 2012 → 2013, com um indicador consistente | −3,00% |
| **2013 → 2014, como o arquivo a mostra** | **−6,70%** |
| 2013 → 2014, com o mesmo indicador dos dois lados | −3,86% |

Cerca de metade da descida aparente vem da mudança de fonte. Numa base
consistente, 2014 continua a tendência anterior e não há nada de especial.

A análise de quebras estruturais assinala este ponto como quebra na série
padronizada de Portugal. **Não é um acontecimento epidemiológico.** Quando uma
quebra detectada coincide com uma mudança de base, a aplicação passou a dizê-lo
no texto da análise.

Ao nível municipal o efeito é maior: 126 dos 305 municípios movem-se mais de 2%
nessa passagem, contra 9 a 19 nos anos vizinhos. Corvo −13,8%, Faro +7,2%,
Lisboa +6,3%.

A mudança de base **não é removível**: o indicador mais recente só começa em
2011, e não existe série municipal anterior na base nova. Pode ser deslocada,
não eliminada.

Nas taxas brutas o efeito é desprezável (0,2%). Se a sua análise atravessa
2013/2014 e precisa de uma taxa comparável, a taxa bruta é a opção segura.

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

- **ULS (Unidade Local de Saúde):** a unidade de organização do SNS a que corresponde a população de um conjunto de municípios. Na aplicação é somada a partir dos seus municípios.
- **ARS:** as cinco regiões de saúde (Norte, Centro, Lisboa e Vale do Tejo, Alentejo, Algarve) que agrupam as ULS. Não coincidem com as regiões NUTS com o mesmo nome.
- **NUTS:** a nomenclatura estatística das regiões. A aplicação usa dois níveis: NUTS I (`Continente`, Açores, Madeira) e NUTS II (as regiões).
- **Definição das regiões (NUTS 2013 / NUTS 2024):** as duas versões da nomenclatura que a aplicação oferece, no controlo do topo da página. Agrupam os mesmos 308 municípios de formas diferentes.
- **Área Metropolitana de Lisboa:** a região de Lisboa em NUTS 2013. Em NUTS 2024 está dividida em `Grande Lisboa` e `Península de Setúbal`.
- **Oeste e Vale do Tejo:** região criada em NUTS 2024 com o Oeste, o Médio Tejo e a Lezíria do Tejo. Apesar do nome, não inclui Lisboa.
- **Agregação por municípios:** as regiões são sempre somadas a partir dos seus municípios, com a mesma lista aplicada a todos os anos. Mantém a série contínua, mas os totais não coincidem exactamente com os publicados pelo INE.
- **Revisão da população (2021):** o INE publica duas estimativas que não coincidem; a aplicação usa a revista a partir de 2021. Todas as taxas desde esse ano mudaram, e há um degrau de cerca de 1,7% em 2020/2021.
