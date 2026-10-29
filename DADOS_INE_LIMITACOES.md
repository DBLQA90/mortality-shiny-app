# Dados do INE: incoerências, limitações e abordagens

Este documento reúne os problemas encontrados nos dados que a aplicação importa do INE, ou na forma como os líamos. Para cada um diz:

- o que é;
- se é **definitivo** (está na fonte e não vai mudar) ou **temporário** (desaparece quando o INE publicar, ou quando nós corrigirmos);
- que dimensão tem;
- que abordagens eram possíveis;
- qual é a usada por omissão, e porquê.

A última validação completa foi feita a 22 de outubro de 2026. As escolhas metodológicas gerais estão em [ESCOLHAS_METODOLOGICAS.md](ESCOLHAS_METODOLOGICAS.md) e os pormenores técnicos em [METHODOLOGY.md](METHODOLOGY.md).

## Resumo

São 32 pontos: 26 definitivos, 4 temporários e 2 que são escolhas de método (denominadores). Quatro pesam de facto nos números, e todos têm tratamento por omissão:

- **A idade dos óbitos nos municípios é incompleta,** sobretudo em 2014 (faltam 31% das idades no cancro do pulmão). Os óbitos em falta são distribuídos pelas idades, e as regiões usam as linhas do INE, que são completas.
- **A população tem duas quebras de série,** em 2013/2014 (mudança de indicador) e em 2020/2021 (revisão, +1,7% a +5,3%). Não se podem remover: ficam onde fazem menos estragos, e a aplicação avisa.
- **O INE publica errados os óbitos de 2011 de cinco municípios do Alto Alentejo** (139 em vez de 424). Usa-se o total por freguesia.
- **O INE oculta valores** para proteger empregadores ou pessoas: totais de óbitos, trabalhadores e ganhos. Recuperam-se quando os totais o permitem (`≈`); quando não, a área fica sem valor.

A validação também encontrou erros nossos, todos corrigidos:

- brancos do INE guardados como zero, nos totais de óbitos e nos indicadores de planeamento;
- municípios em falta num download de 2022–2023;
- leituras pelo nome (Lisboa, Calheta, Lagoa).

O único problema temporário com peso é o ano mais recente sem causas de morte (hoje 2025), que é reconstruído e marcado `p`.

O botão **Dados: corrigidos / originais do INE**, no cabeçalho, permite ver os números com ou sem as correções da aplicação (ver o fim deste documento).

## Como foi validado

Os dados são lidos tal como o INE os publica, sem as correções da aplicação, e cada conjunto é comparado com outro que meça a mesma coisa. A auditoria está em `tools/source_audit.R`, corre em cerca de 10 segundos e pode repetir-se a cada atualização.

| Verificação | O que se compara | Resultado (1991–2025) |
|---|---|---|
| Municípios contra Portugal | Soma dos 308 municípios contra a linha nacional do INE, em óbitos, população, nascimentos e óbitos infantis | População e óbitos infantis exatos em todos os anos; óbitos até 0,11% abaixo, exceto em 2011 (0,30%) |
| Sexos | Ambos os sexos contra homens + mulheres, em cada município | Coerente em todos os municípios e anos, exceto 5 municípios em 2011 |
| Idade contra total | Soma dos óbitos por idade contra o total de todas as idades | Faltam idades em vários anos (ver Óbitos) |
| Freguesias contra municípios | Óbitos por freguesia somados por município contra o total municipal do INE | 8 904 de 8 932 municípios-ano exatos (1996–2024) |
| Sub-regiões contra municípios | Linha NUTS III do INE contra a soma dos seus municípios | Exata, salvo diferenças de 3 a 15 óbitos na Área Metropolitana de Lisboa e o caso de 2011 |
| Semanal contra anual | Soma das semanas contra os óbitos por freguesia do ano civil | Entre −0,2% e −0,7% (2021–2025), pelas semanas ISO |
| Fontes sobrepostas | Indicadores do INE que publicam os mesmos anos (população 2011–2013 e 2021–2023; nascimentos 2022) | Ver População e Nascimentos |
| Valor publicado | Pedido direto à API do INE nos casos suspeitos | Confirmado o erro de 2011 e o estado «Dado confidencial» dos totais em branco |

Há ainda duas verificações que correm a cada atualização:

- **Regressão:** 119 valores validados contra indicadores do INE por município (taxas brutas de natalidade e mortalidade, índice de envelhecimento, mortalidade infantil quinquenal, entre outros).
- **Invariantes:** 342 cálculos em várias áreas, anos e opções, à procura do que nunca deveria acontecer.

## Óbitos

Os óbitos têm os problemas com mais peso: a repartição por idade incompleta nos municípios e um erro de publicação do INE em 2011.

### 1. Repartição por idade incompleta nos municípios

- **O que é.** O INE publica todos os óbitos no total de cada município, mas parte deles não aparece nos grupos de idade: óbitos sem idade registada e células pequenas ocultadas.
- **Tipo e dimensão: definitivo.** Em todas as causas faltam 0,5–1% das idades em 1991–2001, quase nada em 2002–2012, 0,7% em 2013, **3,9% em 2014** (123 municípios perdem mais de 10%), 0,5% em 2015 e 0,4% em 2024. Nas causas específicas é pior: no cancro do pulmão faltam 1–4% na maioria dos anos e **31% em 2014**. Nenhum outro indicador do INE resolve 2014.
- **Abordagens:**
  - usar só os óbitos com idade, o que subestima as taxas;
  - distribuir os que faltam pelo perfil de idade do próprio município, o que desacerta os totais nacionais por idade;
  - distribuir pelo perfil nacional dos que faltam, que é enviesado para os jovens (dava 19 óbitos de menores de 5 anos em Alvito, 1997);
  - **ajustar às duas margens conhecidas** (quantos faltam em cada município e em que idades faltam no país).
- **Por omissão: ajustar às duas margens,** em todos os separadores. É a única abordagem que respeita tudo o que se sabe. Acima de 25% de óbitos sem idade, o separador de planeamento não apresenta taxas por idade.

### 2. Regiões: linha do INE ou soma dos municípios

- **O que é.** Pelo problema anterior, somar municípios perde óbitos por causa. As linhas regionais do INE são completas: as idades somam o total em todos os anos, 2014 incluído.
- **Tipo e dimensão: definitivo.** No cancro do pulmão em 2014, a soma dos municípios perde 31% no Norte e 84% na Madeira.
- **Abordagens:** linhas regionais do INE, ou soma dos municípios.
- **Por omissão: linhas regionais do INE,** com a soma dos municípios disponível no cabeçalho. As regiões redesenhadas em 2024 compõem-se a partir das sub-regiões, e cada composição bate certo com 2022, o único ano publicado nas duas versões.

### 3. Totais confidenciais

- **O que é.** O INE deixa em branco («Dado confidencial») o total de todas as idades de alguns municípios, embora publique as idades.
- **Tipo e dimensão: definitivo.** Em 1991–2024: 8 totais de todas as causas (por exemplo Vimioso 2024, Miranda do Douro e Alfândega da Fé 2015) e 2 787 totais por causa.
- **Abordagens:** tratar como 0 (errado), tratar como em falta, ou usar a soma das idades.
- **Por omissão: soma das idades,** porque um total nunca pode ser menor do que as suas partes. Até 22 de outubro de 2026, os nossos ficheiros guardavam estes brancos como 0, por um erro da nossa extração (o pacote de leitura do INE devolvia-os bem). Os ficheiros foram descarregados de novo; os resultados não mudaram, porque a aplicação já substituía esses zeros.

### 4. Erro do INE em 2011: cinco municípios do Alto Alentejo

- **O que é.** A linha de todas as causas, ambos os sexos, de Nisa, Gavião, Fronteira, Monforte e Marvão está errada tal como o INE a publica (indicador 0008206): 139 óbitos no conjunto, em vez de 424.

  | Município | Publicado pelo INE | Óbitos por freguesia |
  |---|---:|---:|
  | Nisa | 24 | 165 |
  | Gavião | 34 | 91 |
  | Fronteira | 1 | 60 |
  | Monforte | 31 | 49 |
  | Marvão | 49 | 59 |

- **Tipo e dimensão: definitivo,** até o INE o corrigir. Confirmado na API do INE (atualização de 11/10/2025). Três fontes independentes dão 424:
  - os óbitos por freguesia;
  - a linha do INE para o Alto Alentejo (1 888), que nos outros anos é exatamente a soma dos municípios;
  - os próprios sexos, que são incoerentes: em Gavião há 91 homens para 34 óbitos no total.

  Não se trata de óbitos de município desconhecido: todos têm freguesia. As linhas por causa parecem certas.
- **Abordagens:** usar os valores publicados, ou usar o total por freguesia com a repartição por sexo e idade de 2010 e 2012.
- **Por omissão: o total por freguesia.** Os ficheiros guardam os valores do INE tal como publicados, e a correção é aplicada à leitura a partir de `data/source_corrections.rds` (construído por `tools/build_source_corrections.R`). As regiões e a ULS Alto Alentejo usam as linhas regionais e nunca foram afetadas.

### 5. Óbitos de residência desconhecida

- **O que é.** A linha nacional do INE inclui residentes sem município conhecido. A soma dos municípios não os tem.
- **Tipo e dimensão: definitivo.** Entre 0% e 0,11% dos óbitos, consoante o ano.
- **Abordagens:** Portugal = total do INE, ou Portugal = soma dos municípios.
- **Por omissão: as duas, à escolha no separador de planeamento.** Para comparar um local com o país, a soma dos municípios compara igual com igual.

### 6. Pequenas diferenças sem explicação

- **O que é:**
  - a linha do INE para a Área Metropolitana de Lisboa fica 3 a 15 óbitos abaixo da soma dos seus municípios em 12 anos entre 2006 e 2022;
  - no Algarve, em 2010, fica 3 abaixo;
  - em Ponte de Lima 2019, o total municipal tem 469 óbitos e as freguesias 472.
- **Tipo e dimensão: definitivo e desprezável** (menos de 0,06%).
- **Por omissão: valores do INE sem correção.**

### 7. Linhas regionais com nomes repetidos

- **O que é.** No INE, Algarve, Grande Lisboa, Península de Setúbal, Açores e Madeira têm o mesmo nome em dois ou três níveis NUTS. Lidas pelo nome, essas linhas contavam os óbitos duas ou três vezes.
- **Tipo: definitivo na fonte, resolvido.** A aplicação identifica as linhas pelo código geográfico, nunca pelo nome.

### 8. O último ano sem causas de morte

- **O que é.** Os óbitos por causa e idade chegam cerca de 15 meses depois do fim do ano. Hoje falta 2025.
- **Tipo: temporário.**
- **Abordagens:** esperar pelo INE, ou reconstruir o ano a partir dos óbitos semanais (que dão as idades) e dos óbitos por freguesia (que dão o total).
- **Por omissão: reconstruir,** só para Portugal, o Continente e as regiões NUTS II inteiras, em todas as causas e ambos os sexos, com a marca `p`. Testado em 2022–2024: diferenças até 0,13% nos óbitos, 0,03% na taxa padronizada e 0,3% na prematura. As semanas ISO, somadas por ano, ficam 0,2–0,7% abaixo do ano civil, e por isso o nível vem das freguesias.

### 9. Lista de causas

- **O que é.** A lista sucinta do INE (0008206, 0013166) junta o cancro da laringe ao da traqueia, brônquios e pulmão (a laringe é cerca de 8%). O indicador 0008461 separa-os, e detalha os tumores malignos em 60 localizações, mas só cobre 2002–2022.
- **Tipo e dimensão: definitivo.**
- **Por omissão: a lista sucinta,** que cobre toda a série. As 60 localizações do 0008461 estão disponíveis como causas adicionais, só nos anos 2002–2022.
- **No 0008461, a soma das localizações fica abaixo do total em 2002–2018,** 7 a 34 óbitos por ano no Continente (0,03–0,13%); a partir de 2019 coincide. O total coincide com o do 0008206. A diferença é demasiado pequena para mudar uma taxa e não é corrigida; foi perguntada ao INE.

### 10. Quebras de codificação das causas de morte

- **O que é.** Em algumas causas, a taxa muda de nível de um ano para o outro sem que a mortalidade o explique: muda a forma como os óbitos são codificados.
  - **Transtornos mentais e comportamentais:** a taxa padronizada passa de cerca de 1–2 por 100 000 até 2012 para cerca de 21 em 2013, e continua a subir (Continente: 1,7 em 2012, 20,8 em 2013, 46,3 em 2023). Coincide com a mudança na codificação das demências, com o início do SICO.
  - **Doenças da pele e do tecido celular subcutâneo:** a taxa duplica de 2018 para 2019 (Região Norte: 1,7 para 3,4).
- **Tipo e dimensão: definitivo.** As séries não são comparáveis antes e depois da quebra.
- **Abordagens:** ignorar a quebra (as projeções extrapolam o salto: no Norte, mulheres, transtornos mentais, 177,8 em 2030 contra 53,7 em 2023), usar só os anos depois da quebra, ou não projetar.
- **Por omissão:** a aplicação ainda não deteta estas quebras; a série observada mostra-as tal como estão. Nas projeções de 2030 da tabela «82_MRT» (outubro de 2026): transtornos mentais ajustados só a 2013–2024; doenças da pele sem projeção, porque depois de 2019 há só 6 anos. Detetar estas quebras na Previsão Guiada e propor o ajuste depois delas é uma melhoria prevista.

## População

A população é internamente coerente: em todos os anos de 1991 a 2025, os 308 municípios somam exatamente Portugal, e homens mais mulheres somam o total. Os problemas estão na passagem entre indicadores do INE, que não concordam nos anos que têm em comum.

### 1. Mudança de indicador em 2013/2014

- **O que é.** Até 2013, a população vem do indicador 0003182 (NUTS 2002); de 2014 a 2020, do 0008273. Os dois distribuem a população pelas idades de forma diferente.
- **Tipo e dimensão: definitivo.** O total quase não muda (0,16%), mas a taxa padronizada de Portugal fica 2,5–3% diferente. Cerca de metade da descida aparente de 2013 para 2014 vem só da mudança. Em 2014, 66 municípios mudam mais de 3% de um ano para o outro, contra 1 a 10 nos anos de 1992 a 2013.
- **Abordagens:** mudar de indicador noutro ano, entre 2011 e 2013 (o salto ficaria parecido), ou manter e avisar. Não existe série municipal única de 1991 a 2020.
- **Por omissão: mudar em 2014 e avisar** sempre que uma série, tendência ou projeção atravessa esse ano. É o ponto onde o indicador antigo acaba. Com isto, 2011–2013 ficam 0,16% abaixo do valor atual do INE.

### 2. Revisão da população a partir de 2021

- **O que é.** O INE publica duas estimativas para os mesmos anos: a antiga (0008273, até 2023) e a revista com base administrativa (0012918, 2021–2025). A revista é mais alta, e a diferença cresce de ano para ano.
- **Tipo e dimensão: definitivo.** A diferença é de +1,7% em 2021, +3,9% em 2022 e +5,3% em 2023. As mudanças de mais de 3% num ano tocam 92 municípios em 2021, 79 em 2022 e 46 em 2023.
- **Abordagens:**
  - usar a antiga até 2023 e a revista depois: dá um degrau de 5,3% em 2024, que parece uma queda da mortalidade;
  - **usar a revista em todo o seu período (desde 2021):** o degrau fica em 2020/2021, com 1,7%.
- **Por omissão: a revista desde 2021, com aviso.** Todas as taxas desde 2021 desceram (a mortalidade bruta de Portugal em 2023 é 1 056 por 100 000, e era 1 112 com a série antiga). Antes de 2021 não há revisão: a série nacional longa do INE coincide com a aplicação, pessoa a pessoa.

### 3. «Lisboa», «Calheta» e «Lagoa»

- **O que é.** Na população até 2013, «Lisboa» era ao mesmo tempo a região e o município, e em alguns indicadores há duas Calhetas e duas Lagoas (uma nos Açores, outra na Madeira ou no Algarve). Lidas pelo nome, juntavam-se. A taxa de Lisboa ficava 6 vezes abaixo do real até 2013.
- **Tipo: definitivo na fonte, resolvido.** Os dados são pedidos pelo código geográfico.

### 4. Denominador

- **O que é.** As estimativas do INE referem-se a 31 de dezembro, mas as taxas do INE dividem pela média de dois fins de ano.
- **Tipo:** não é um problema, é uma escolha.
- **Por omissão: a média,** com 31 de dezembro disponível no cabeçalho. Com a média, a aplicação reproduz as taxas brutas do INE nos 308 municípios.

### 5. Grupo aberto dos 85 e mais anos

- **O que é.** A população por município para no grupo «85 e mais anos».
- **Tipo: definitivo.** Limita a esperança de vida: a tábua abreviada fica cerca de 0,8 anos acima da do INE, que usa idades simples, e de acordo com a do Eurostat.

## Nascimentos e óbitos infantis

Hoje, os nascimentos e os óbitos infantis batem certo com Portugal. Os problemas foram de leitura, já resolvidos, e de cobertura antes de 2011, que é definitiva. A soma dos municípios fica 0 a 18 nascimentos abaixo de Portugal (residência desconhecida, até 2009) e é exata desde 2010 (salvo 1 nascimento em 2019). Nos óbitos infantis é exata em 2011–2025, com 1 óbito de diferença em 2024.

### 1. Nascimentos lidos pelo nome

- **O que é.** No indicador de nascimentos de 1995–2013 (0000003), «Lisboa» designa a região e o município. Lidos pelo nome, o município recebia 37 208 nascimentos em 2001, quando teve 5 604.
- **Tipo: definitivo na fonte, resolvido.** Leitura pelo código do município; um ano só é aceite se os 308 municípios tiverem correspondência.

### 2. Óbitos infantis antes de 2011

- **O que é.** Até 2010, a única fonte municipal é o grupo «menos de 1 ano» dos óbitos por causa, que sofre da mesma repartição por idade incompleta (em 2014: 112 óbitos nos municípios, contra 236 no país). Desde 2011, o INE publica os óbitos infantis como tema próprio (0008181 e 0012541), e esses são completos.
- **Tipo e dimensão: definitivo.** A soma dos municípios atinge cerca de 85% de Portugal em 1995–2001 e 99–100% em 2002–2010.
- **Abordagens:** não mostrar esses anos, ou mostrá-los assinalados.
- **Por omissão: mostrar desde 2011 a fonte completa** e, antes disso, assinalar com `†` os anos abaixo de 97%. Portugal usa a sua própria linha e não é afetado.

### 3. Municípios com poucos nascimentos

- **O que é.** Muitos municípios têm menos de 10 nascimentos por ano. Um óbito muda a taxa em centenas por 1 000.
- **Tipo: definitivo** (é a realidade, não um defeito da fonte).
- **Abordagens:** ocultar os valores com poucos nascimentos, ou mostrá-los assinalados.
- **Por omissão: mostrar sempre,** com `*` abaixo de 1 000 nascimentos no período e o número de óbitos infantis como indicador próprio. Nada é ocultado.

### 4. Fetos mortos

- **O que é.** O INE publica os óbitos perinatais por município, mas não os fetos mortos isolados.
- **Tipo: definitivo.**
- **Por omissão:** fetos mortos = óbitos perinatais menos óbitos com menos de 7 dias. Por isso, a mortalidade fetal tardia e a perinatal começam no triénio 2011–2013.

### 5. 2025

- **O que é.** Os óbitos infantis de 2025 só existem em todas as causas e para ambos os sexos (0012540).
- **Tipo: temporário.**
- **Por omissão:** usar 2025 onde o indicador o permite (mortalidade infantil de todas as causas).

## Geografia

Os limites dos 308 municípios não mudaram desde 1991. O que muda são os agrupamentos acima deles, e os códigos e nomes com que o INE os identifica.

### 1. Duas versões das NUTS

- **O que é.** Em 2024, o INE redesenhou as regiões. A Lezíria do Tejo saiu do Alentejo, e Sertã e Vila de Rei passaram para a Beira Baixa. Lidas sem cuidado, as linhas regionais comparam dois Alentejos: 11 327 óbitos em 2022 numa versão e 7 898 na outra.
- **Tipo: definitivo.**
- **Abordagens:** usar a definição em vigor em cada ano, ou **uma lista fixa de municípios para todos os anos**.
- **Por omissão: lista fixa, NUTS 2024,** com NUTS 2013 disponível no cabeçalho. A série fica contínua. Uma região é «este território», e não «o que se publicou com este nome em cada ano».

### 2. Códigos que mudam

- **O que é.** O código de cada município inclui as suas NUTS. Por isso, só 132 dos 308 mantêm o código entre versões.
- **Tipo: definitivo, resolvido.** A aplicação liga os municípios pelos últimos quatro dígitos do código (DICO), que nunca mudam. Quando não há correspondência, dá erro em vez de adivinhar.

### 3. Municípios criados em 1998

- **O que é.** Odivelas, Trofa e Vizela foram criados em 1998, a partir de Loures, Santo Tirso e Guimarães. O INE recalculou a população para os limites atuais, mas os nascimentos e os óbitos até 1998 ficaram no município de origem. Sem cuidado, a Trofa tinha uma esperança de vida de 92,4 anos.
- **Tipo: definitivo.**
- **Por omissão:** até 1998, qualquer área que tenha um dos municípios do par sem o outro fica sem valor. O par junto (por exemplo, Loures com Odivelas) tem valor.

### 4. ULS que partilham um município

- **O que é.** Lisboa, Loures e Porto estão divididos entre ULS por freguesia. Gondomar também, por uma só freguesia, a Lomba (0,8% dos residentes).
- **Tipo: definitivo.**
- **Abordagens:** dar o município inteiro a cada ULS, ou repartir pelas freguesias.
- **Por omissão: as duas, à escolha no separador de planeamento.** Com o município inteiro, Gondomar fica todo na ULS Santo António. A repartição usa os nascimentos e os óbitos registados por freguesia, que o INE publica todos os anos desde 2014. A população reparte-se pelas quotas das freguesias nos Censos de 2021, porque não há estimativas anuais por freguesia.

### 5. Reorganização das freguesias em 2013

- **O que é.** Até 2013, os óbitos por freguesia vêm com as freguesias antigas, e as «uniões de freguesias» só existem depois. Nesses anos, as uniões aparecem vazias.
- **Tipo: definitivo.** Não é falta de dados: as freguesias antigas somam o município.
- **Por omissão:** antes de 2014, a repartição das ULS usa as quotas do Censo.

### 6. Municípios em falta num download (erro nosso)

- **O que é.** Os óbitos de 2022 e 2023 tinham sido descarregados município a município, com uma lista desatualizada, e faltavam as duas Calhetas e a Lagoa dos Açores. Os Açores ficavam 7% abaixo.
- **Tipo: temporário, resolvido em agosto de 2026.** Desde então, cada ano só é aceite se os municípios somarem a linha nacional.

## Outros indicadores

Nos indicadores socioeconómicos, o problema típico é o INE ocultar o valor de um município para proteger um empregador ou uma pessoa. Uma ocultação não é zero.

### 1. Emprego e ganhos ocultados (Quadros de Pessoal)

- **O que é.** O INE oculta os trabalhadores ou o ganho de alguns municípios: os 19 municípios dos Açores em 2013–2014 (trabalhadores) e 2011–2014 (ganhos), São João da Madeira em 2015, e outros. Nos setores de atividade, oculta dois de três quando um só revelaria uma empresa (6 a 49 municípios por ano).
- **Tipo: definitivo.**
- **Abordagens:** área sem valor, ou recuperar o valor.
- **Por omissão: recuperar,** com a marca `≈`. Os totais do Continente e de Portugal ainda incluem o valor ocultado, e por isso o valor de um município é o total menos os municípios publicados. Nos anos sem ocultação, este método reproduz o total exato nos trabalhadores e a menos de 0,005% nos ganhos. Os setores ocultados repartem-se como na restante NUTS III.

### 2. Resíduos e pensões em falta

- **O que é:**
  - faltam resíduos de muitos municípios das ilhas em 1995–2006;
  - Loures e Odivelas têm um serviço comum (SIMAR), registado todo em Loures;
  - faltam pensionistas de 4 municípios em 2005;
  - as pensões mudam de série em 2017, com cerca de 5,5% menos pensionistas.
- **Tipo: definitivo.** Aqui os totais do INE também deixam de fora os municípios em falta, e não dá para recuperar.
- **Por omissão:** a área que contém um município sem valor fica sem valor. Os resíduos por habitante só existem para áreas com Loures e Odivelas juntos. A quebra de 2017 está marcada no gráfico.

### 3. Denominadores dos indicadores do INE

- **O que é.** As taxas de beneficiários do RSI e de pensionistas do INE dividem pela população com 15 ou mais anos.
- **Por omissão: 15 e mais anos,** como o INE. A aplicação fica a cerca de 1% dos valores do INE por município em 2024.

### 4. Censos

- **O que é.** A escolaridade por idade e município só existe nos Censos de 2011 e 2021, e os censos usam códigos de município com outro formato.
- **Tipo: definitivo.**
- **Por omissão:** com idade mínima, 1991 e 2001 ficam sem valor.

### 5. Cuidados de saúde primários (Portal da Transparência do SNS)

- **O que é.** Muitos indicadores acumulam ao longo do ano e recomeçam: os rastreios e o pé diabético por ano civil, a tensão arterial e a HbA1c por semestre. Os dados são por ULS desde 2024 e por ACES antes.
- **Tipo: definitivo.**
- **Por omissão:** comparar dezembro (ou junho, nos semestrais). Um valor de julho não é comparável.

### 6. Brancos guardados como zero (erro nosso)

- **O que é.** O script que descarrega estes indicadores somava as componentes de cada município de forma que um valor em branco do INE ficava 0. Por exemplo, os trabalhadores dos Açores em 2013 estavam guardados como 0.
- **Tipo: temporário, resolvido a 22 de outubro de 2026.** O script foi corrigido e os 19 conjuntos foram descarregados de novo; os brancos ficam agora em falta.
- **O que se encontrou.** Só há brancos onde já se conheciam:
  - trabalhadores por setor (685 células em 12 anos) e ganhos (os 19 municípios dos Açores em 2011–2014);
  - resíduos (481 células em 27 anos, sobretudo ilhas em 1995–2006);
  - pensão média (4 municípios em 2005);
  - censo de 1991 e poder de compra nos anos anteriores a 1998 (Odivelas, Trofa e Vizela, que ainda não existiam).

  Nos nascimentos por peso, por gestação e por idade da mãe, nos óbitos perinatais, nos óbitos infantis por idade, nos pensionistas e no RSI não há nenhum branco.
- **Efeito: nenhum nos resultados.** Todos estes casos já eram tratados pela aplicação.

## Botão «Dados: corrigidos / originais do INE»

No cabeçalho, ao lado de «Definição das regiões» e de «Óbitos das regiões», o controlo **Dados** tem duas opções: **Corrigidos (recomendado)**, a predefinição, e **Originais do INE**. A segunda desliga tudo o que a aplicação acrescenta ou substitui nos dados do INE.

| O que a aplicação faz | Com «corrigidos» | Com «originais do INE» |
|---|---|---|
| Erros do INE (2011, Alto Alentejo) | Total por freguesia | Valores publicados |
| Óbitos sem idade | Distribuídos pelas idades | Só os óbitos com idade; os totais e as taxas brutas continuam completos |
| Totais confidenciais | Soma das idades | Em falta |
| Linhas regionais sem os óbitos sem idade | Completadas até ao total dos municípios | Como publicadas |
| Emprego e ganhos ocultados | Recuperados dos totais (`≈`) | Em falta |
| Setores de emprego ocultados | Estimados pela NUTS III (`≈`) | Em falta |
| Último ano sem causas (2025) | Reconstruído (`p`) | Só onde o INE publica |

O botão não desliga, porque não são correções ao INE mas escolhas de método, ou correções de erros nossos:

- a definição das regiões e a origem dos óbitos regionais, que têm controlo próprio;
- o denominador (média ou 31 de dezembro), que também tem controlo próprio;
- a série de população revista desde 2021: a série antiga também é «original», e a escolha entre as duas é de método;
- a leitura por código geográfico (Lisboa, Calheta, Lagoa, linhas regionais repetidas): ler pelo nome dá dados que o INE nunca publicou;
- os municípios criados em 1998 e os conjuntos completos: mostrar 0 onde o INE não publicou nada seria inventar um valor.

## Questões enviadas ao INE

As incoerências que só o INE pode esclarecer (o erro de 2011, as diferenças da Área Metropolitana de Lisboa e de Ponte de Lima, a confidencialidade dos totais, a repartição por idade de 2014 e as séries de população sobrepostas) estão num relatório à parte, preparado em outubro de 2026 para envio ao INE.
