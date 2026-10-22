# Escolhas metodológicas, explicadas

Para calcular um indicador, é preciso tomar decisões que o número final não mostra. De onde vêm os óbitos de uma região? Que população se usa no denominador? O que fazer quando o INE não publica um valor? Como se projecta o futuro?

Este documento junta essas decisões, em linguagem simples. Para cada uma diz o que a aplicação faz, porquê, e o que isso muda na leitura dos números. Não é preciso lê-lo de seguida: cada escolha vale por si.

Os pormenores técnicos, com as fórmulas e as verificações, estão em [METHODOLOGY.md](METHODOLOGY.md) (em inglês). O uso de cada separador está no [Manual do Utilizador](USER_MANUAL.md).

**Índice**

1. [De onde vêm os dados](#1-de-onde-vêm-os-dados)
2. [Geografia: regiões, ULS e Portugal](#2-geografia-regiões-uls-e-portugal)
3. [Como se calculam as taxas](#3-como-se-calculam-as-taxas)
4. [Dados incompletos e erros nas fontes](#4-dados-incompletos-e-erros-nas-fontes)
5. [Comparar áreas e ler a incerteza](#5-comparar-áreas-e-ler-a-incerteza)
6. [Tendências e projecções](#6-tendências-e-projecções)
7. [Como se verifica que os números estão certos](#7-como-se-verifica-que-os-números-estão-certos)

---

## 1. De onde vêm os dados

### Dados oficiais, guardados com data

**O que fazemos.** Os dados vêm do INE: óbitos por causa, idade, sexo e município de residência; estimativas de população; nados-vivos; óbitos por freguesia e por semana. Os dados dos cuidados de saúde primários vêm do Portal da Transparência do SNS. A aplicação guarda uma cópia local e actualiza-a todas as semanas. Quando o INE revê um valor, a versão anterior é arquivada, e cada importação fica registada com a data.

**Porquê.** O INE revê dados sem aviso. Com o registo, sabe-se sempre com que versão foi calculado um número, e é possível refazer um cálculo antigo.

**O que significa para si.** Um valor pode mudar de uma semana para a outra se o INE o tiver revisto. A data dos dados aparece nas exportações.

### Local de residência, não de óbito

**O que fazemos.** Os óbitos contam no município onde a pessoa residia, não onde morreu.

**Porquê.** É a convenção do INE e a única que permite calcular taxas: o denominador é a população residente.

**O que significa para si.** Um município com um hospital central não aparece com mais mortalidade por receber doentes de fora.

---

## 2. Geografia: regiões, ULS e Portugal

### As regiões são sempre os mesmos municípios

**O que fazemos.** Uma região (NUTS II ou III) é definida pela lista dos seus municípios, e essa lista é a mesma em todos os anos da série. No topo da página escolhe-se a definição: NUTS 2024 (predefinição) ou NUTS 2013.

**Porquê.** Em 2024 o INE redesenhou várias regiões. A Lezíria do Tejo, por exemplo, passou do Alentejo para o Oeste e Vale do Tejo. Se cada ano usasse a definição em vigor nesse ano, a série compararia dois Alentejos diferentes. Com a lista fixa, a série é contínua.

**O que significa para si.** Atenção aos nomes iguais: o `Centro` e o `Alentejo` não são o mesmo território em NUTS 2013 e em NUTS 2024. A definição activa está sempre à vista no topo da página.

### Óbitos das regiões: a linha do próprio INE

**O que fazemos.** Por omissão, os óbitos de uma região vêm da linha que o INE publica para essa região. Quando a região não tem linha própria num ano, compõe-se a partir das sub-regiões. A população é sempre a soma dos municípios.

**Porquê.** Nos municípios, o INE publica os totais completos, mas a repartição por idade fica incompleta, sobretudo onde os números são pequenos. Como as taxas padronizadas precisam da idade, somar municípios perderia óbitos: no cancro do pulmão, em 2014, a perda chegou a 30–84% em todas as regiões. As linhas regionais do INE não têm esse problema.

**O que significa para si.** Os valores regionais são fiáveis por causa e por idade. Um município isolado, ou um conjunto de municípios escolhido por si, não tem essa protecção, e a aplicação avisa quando há risco de subestimar.

### Região-plano de Lisboa e Vale do Tejo

**O que fazemos.** A região-plano LVT é a soma de três NUTS II: Grande Lisboa, Península de Setúbal e Oeste e Vale do Tejo. Somam-se os óbitos e as populações e só depois se calcula a taxa.

**Porquê.** A taxa de uma área maior não é a média das taxas das suas partes. Com a média, uma parte pequena pesaria tanto como uma grande.

**O que significa para si.** Vale para qualquer agregação feita na aplicação: seleccionar várias áreas soma-as primeiro e calcula depois.

### ULS que partilham um município

**O que fazemos.** Lisboa, Loures e Porto estão divididos entre ULS ao nível da freguesia. Para estas ULS, o separador de planeamento dá duas leituras à escolha:

- **município inteiro:** cada ULS leva o município todo;
- **ponderação por freguesias:** cada ULS leva só as suas freguesias.

Nos outros separadores conta o município inteiro. Gondomar fica todo na ULS Santo António, embora a freguesia da Lomba (0,8% dos residentes) seja da ULS Vila Nova de Gaia/Espinho. A Lomba só passa para Gaia/Espinho com a ponderação por freguesias.

**Porquê.** A leitura certa depende da pergunta. Uma ULS não pode recusar residentes do município que serve, e para isso conta o município inteiro. Para comparar ULS entre si, ou somá-las sem contar ninguém duas vezes, contam as freguesias. Dar Gondomar inteiro a Gaia/Espinho por causa de uma freguesia somar-lhe-ia 164 mil residentes.

**O que significa para si.** Com município inteiro, estas ULS sobrepõem-se: somá-las todas conta Lisboa, Loures e Porto mais de uma vez. Ao escolher uma delas, a aplicação diz como está a ser lida.

### Portugal: total do INE ou soma dos municípios

**O que fazemos.** No separador de planeamento, Portugal pode ser o total publicado pelo INE ou a soma dos 308 municípios.

**Porquê.** O total do INE inclui residentes cujo município é desconhecido, entre 0,3% e 0,9% dos óbitos. Para perguntar se um local se afasta do país, compara-se igual com igual, e por isso usa-se a soma dos municípios.

**O que significa para si.** A diferença é pequena, mas pode decidir um teste de significância no limite.

---

## 3. Como se calculam as taxas

### Denominador: a população média do ano

**O que fazemos.** As taxas dividem os acontecimentos de um ano pela média da população no fim do ano anterior e no fim desse ano. É uma aproximação da população a meio do ano. No topo da página pode escolher-se, em alternativa, a população a 31 de Dezembro.

**Porquê.** É assim que o INE calcula as suas taxas. Com esta escolha, a aplicação reproduz as taxas brutas do INE em todos os municípios.

**O que significa para si.** Com a população a 31 de Dezembro, as taxas saem mais baixas onde a população cresce e mais altas onde diminui. Nas regiões, entre 2019 e 2025, a diferença foi de −0,6% a +3,9%.

### Taxa bruta e taxa padronizada

**O que fazemos.** A taxa bruta divide os óbitos pela população. A taxa padronizada (método directo) calcula quantos óbitos haveria se cada área tivesse a mesma estrutura etária: a população padrão europeia de 2013, em 18 grupos de idade, dos 0–4 aos 85 e mais anos.

**Porquê.** Morre mais gente onde há mais idosos. A taxa bruta mistura duas coisas, o risco de morrer e a idade da população. A padronizada separa-as e permite comparar áreas e anos.

**O que significa para si.** Para comparar áreas, use a taxa padronizada. Por exemplo, em 2025, o Alentejo teve uma taxa bruta 1,49 vezes a de Lisboa e Vale do Tejo; padronizada, 1,17 vezes. A maior parte da diferença vem do envelhecimento, não de maior risco.

### Mortalidade prematura

**O que fazemos.** A mortalidade prematura é a taxa padronizada só com as idades abaixo dos 75 anos.

**Porquê.** Concentra-se nas mortes que mais dependem da prevenção e dos cuidados de saúde.

### Juntar vários anos nas áreas pequenas

**O que fazemos.** Qualquer indicador pode ser calculado em janelas de 3 ou 5 anos. Somam-se os óbitos e a população de todos os anos da janela e só depois se divide. Não se faz a média das taxas anuais.

**Porquê.** Num município pequeno, um ano mau pode duplicar a taxa. Juntar anos dá um valor mais estável.

**O que significa para si.** Ganha-se estabilidade e perde-se pormenor: um valor de 5 anos responde a «qual é o nível habitual?», não a «o que aconteceu neste ano?».

### Mortalidade infantil: mostrar sempre, assinalar quando é frágil

**O que fazemos.** A mortalidade infantil é o número de óbitos com menos de 1 ano por 1.000 nados-vivos. Nunca se esconde um valor. Quando o período tem menos de 1.000 nascimentos, o valor leva um asterisco. O número de óbitos infantis existe também como indicador próprio.

**Porquê.** Abaixo de 1.000 nascimentos, um único óbito muda a taxa em mais de 1 por 1.000, um valor maior do que a diferença entre a maioria das regiões. Esconder o valor tiraria informação. O asterisco avisa que comparar ou ordenar esses valores é ler ruído.

### Esperança de vida

**O que fazemos.** A esperança de vida é calculada com uma tábua de mortalidade abreviada, por grupos de 5 anos (método de Chiang II), com intervalo de confiança.

**Porquê.** O INE usa tábuas por idade simples, com métodos próprios para as idades mais avançadas. Esses métodos não são reprodutíveis com os dados que publica por região. A tábua abreviada é o método de referência em saúde pública, permite calcular qualquer área, incluindo a região-plano LVT e as ULS, e dá o intervalo de confiança.

**O que significa para si.** Os valores ficam cerca de 0,8 anos acima dos do INE, e de acordo com os do Eurostat. Compare valores da aplicação entre si, não com os publicados pelo INE.

---

## 4. Dados incompletos e erros nas fontes

### Óbitos sem idade publicada

**O que fazemos.** Quando a soma por idade de um município fica abaixo do seu total, os óbitos em falta são distribuídos pelas idades. Sabem-se duas coisas: quantos óbitos faltam em cada município, e em que idades faltam no país. A distribuição parte da população por idade do município, com as taxas do país, e é ajustada até respeitar as duas. Isto acontece em todos os separadores. No separador de planeamento, um valor com mais de 2% de óbitos sem idade é assinalado, e com mais de 25% não é apresentado.

**Porquê.** O INE publica todos os óbitos no total, mas nem todos têm a idade registada, e as células pequenas são ocultadas. Sem este passo, Barrancos teria 19 óbitos em 2024 em vez de 32, e o mesmo município daria valores diferentes em separadores diferentes.

**O que significa para si.** Os totais coincidem com os do INE. A repartição por idade de uma pequena parte dos óbitos é estimada.

### Contagens a partir dos totais

**O que fazemos.** Números de óbitos, taxas brutas e proporções vêm do total publicado para cada município e causa, não da soma por idade.

**Porquê.** Esses indicadores não precisam da idade, e o total é completo.

### Valores confidenciais: em falta, nunca zero

**O que fazemos.** Quando o INE marca um valor como confidencial, ele fica guardado como em falta. Se o total de um município for confidencial mas as idades forem públicas, usa-se a soma por idade.

**Porquê.** Um valor confidencial não é zero. Em Outubro de 2026 encontrámos ficheiros de totais antigos que guardavam os confidenciais como 0, por um erro nosso na extracção (o pacote usado para ler o INE devolvia-os bem). O script de extracção já tinha sido corrigido, e os ficheiros foram descarregados de novo. Os valores validados não mudaram, porque a aplicação já substituía esses zeros pela soma por idade.

### Emprego e ganhos ocultados

**O que fazemos.** Quando o INE oculta o emprego ou o ganho de um município, recupera-se o valor a partir dos totais do Continente e de Portugal, que o incluem. Os valores recuperados levam a marca `≈`.

**Porquê.** Sem isto, a região que contém esse município ficaria sem valor, ou pior, somaria o município como zero.

### O último ano antes das causas de morte (marca `p`)

**O que fazemos.** O INE publica os óbitos por causa e idade com cerca de 15 meses de atraso. Entretanto, o ano em falta (hoje 2025) é reconstruído para Portugal, o Continente e as regiões NUTS II inteiras, só para todas as causas e ambos os sexos. A distribuição por idade vem dos óbitos semanais, e o total vem dos óbitos por freguesia, que coincidem com os totais definitivos do INE. Estes valores levam a marca `p` (provisório).

**Porquê.** Sem isto, as regiões ficariam um ano mais atrás. O método foi testado em 2022–2024, comparando com os valores definitivos: as diferenças foram de até 0,13% nos óbitos, 0,03% na taxa padronizada e 0,3% na mortalidade prematura.

**O que significa para si.** Municípios, ULS, causas específicas e cada sexo em separado só têm esse ano quando o INE publicar as causas de morte. O valor provisório é substituído automaticamente nessa altura.

### Mudanças nas fontes que parecem mudanças reais

**O que fazemos.** A aplicação avisa quando uma série, uma tendência ou uma projecção atravessa um destes dois pontos:

- **2013/2014:** muda o indicador de população do INE, que reparte a população pelas idades de outra forma. Em Portugal, cerca de metade da descida aparente da taxa padronizada nesse ano vem só desta mudança.
- **2020/2021:** o INE reviu a população a partir de 2021, com uma nova base administrativa (+1,7% em 2021, a crescer nos anos seguintes). As taxas desde 2021 ficaram mais baixas do que seriam na série antiga.

**Porquê.** Nenhuma das duas mudanças se pode remover: não existe uma série municipal única que cubra os dois lados. O mais honesto é dizer onde estão.

**O que significa para si.** Uma descida em 2014 ou em 2021 pode não ser real. Leia o aviso antes de a interpretar.

### Erros nos dados publicados pelo próprio INE

**O que fazemos.** Os dados são verificados contra fontes independentes. Quando o próprio INE publica um valor que não pode estar certo, a aplicação documenta o problema e diz que áreas são afectadas.

**O caso conhecido: 2011, cinco municípios do Alto Alentejo.** No ficheiro de óbitos por causa, a linha «todas as causas, ambos os sexos» de Nisa, Gavião, Fronteira, Monforte e Marvão tem 139 óbitos no conjunto, em vez de 424. Confirmámos directamente no INE que é assim que está publicado. Há três sinais de que o valor está errado:

- os óbitos por freguesia dão 424;
- a linha do INE para a sub-região Alto Alentejo dá 1.888, e a soma dos seus municípios só bate certo com 424;
- nesses municípios, o número de homens é maior do que o total (Gavião: 91 homens, 34 no total).

**O que fazemos com ele.** Desde 22 de Outubro de 2026, o total destes cinco municípios em 2011 é o dos óbitos por freguesia (424), repartido por sexo e idade como em 2010 e 2012. Em Nisa, onde os homens e as mulheres publicados somam o total certo, mantêm-se. Os ficheiros guardam os valores do INE tal como publicados, e a correcção aplica-se ao ler, a partir de um registo à parte (`data/source_corrections.rds`, construído por `tools/build_source_corrections.R`). As linhas por causa não são corrigidas, porque parecem certas.

**O que significa para si.** As regiões (Alentejo, Alto Alentejo) e a ULS Alto Alentejo usam as linhas regionais do INE, que estão certas, e nunca foram afectadas. Nos cinco municípios, a taxa de 2011 deixa de sair cerca de 85% abaixo do real.

---

## 5. Comparar áreas e ler a incerteza

### Intervalos de confiança

**O que fazemos.** Taxas e proporções levam um intervalo de confiança a 95%. Para taxas de acontecimentos usa-se o intervalo exacto de Poisson; para proporções, o binomial.

**Porquê.** Mesmo com dados completos, o número de óbitos de um ano tem uma parte de acaso. Numa área pequena, essa parte é grande.

**O que significa para si.** Um intervalo largo não é um defeito: diz quanto se pode confiar no valor. Os índices de estrutura da população não têm intervalo, porque o INE não publica o seu erro.

### Diferente de Portugal?

**O que fazemos.** Cada valor com intervalo leva uma marca:

- `▲` quando o intervalo fica todo acima de Portugal;
- `▼` quando fica todo abaixo;
- `=` quando inclui o valor de Portugal.

**Porquê.** É o critério usado no PHE Fingertips (Inglaterra). Uma diferença só conta quando o acaso não a explica.

**O que significa para si.** A marca não diz se a diferença é boa ou má. Uma esperança de vida acima de Portugal é boa; uma mortalidade acima não é.

### Funil

**O que fazemos.** O gráfico em funil mostra todas as áreas de uma vez: o valor de cada uma contra o seu tamanho. Duas linhas marcam o que uma área desse tamanho mostraria só por acaso à volta do valor de Portugal, a 95% e a 99,8%.

**Porquê.** As áreas pequenas variam muito por acaso. O funil mostra quais se afastam mais do que o acaso explica.

**O que significa para si.** Com 308 municípios, cerca de 15 ficam fora da linha de 95% só por acaso. Os pontos fora da linha de 99,8% são os que merecem atenção.

### Suavização das áreas pequenas (opcional)

**O que fazemos.** Na opção «Suavizados (Bayes empírico)», o valor de cada área aproxima-se da média de todas, e tanto mais quanto mais pequena é a área.

**Porquê.** Ordenar municípios pelos valores observados põe nos extremos os mais pequenos, por acaso. Os valores suavizados ordenam melhor.

**O que significa para si.** Use os valores suavizados para ordenar e comparar áreas pequenas. Use os observados nos relatórios e para saber se uma área difere de Portugal.

### Mapas: a cor diz se é melhor ou pior

**O que fazemos.** Nos mapas, cada indicador aparece face a Portugal (Portugal = 100), a vermelho onde a área está pior e a azul onde está melhor. Os indicadores sem melhor nem pior, como a natalidade ou a estrutura etária, usam violeta e amarelo.

**Porquê.** Assim, quatro indicadores em unidades diferentes lêem-se lado a lado. Uma mortalidade alta e uma esperança de vida baixa ficam ambas a vermelho.

---

## 6. Tendências e projecções

### Tendência: variação percentual anual

**O que fazemos.** A tendência mede quanto um indicador muda, em média, por ano, com intervalo de confiança. Pode ter até dois pontos de viragem (método Joinpoint, do National Cancer Institute dos EUA).

**Porquê.** Uma percentagem por ano compara-se entre indicadores e áreas. Os pontos de viragem evitam que uma mudança de direcção desapareça numa média.

**O que significa para si.** Quando o período atravessa 2020/2021, a aplicação avisa: a pandemia e a revisão da população caem nesses anos.

### Projecção: a média de 6 modelos

**O que fazemos.** Por omissão, a Previsão Guiada apresenta a média de seis modelos de séries temporais: ARIMA, ETS, Theta, TBATS, Holt e Holt amortecido. O intervalo junta a incerteza de cada modelo com o desacordo entre eles.

**Porquê.** Até Outubro de 2026 a aplicação escolhia «o melhor modelo», o que errava menos a prever os últimos anos conhecidos. Testámos essa escolha com dados passados: a partir de cada ano entre 2008 e 2020, comparámos a projecção a 4 anos com o que se observou depois, em 27 séries regionais e em 40 municípios e ULS.

| A 4 anos | Regiões: erro mediano | Regiões: o intervalo de 95% contém o valor | Municípios e ULS: erro mediano | Municípios e ULS: contém o valor |
|---|---:|---:|---:|---:|
| O melhor modelo (método anterior) | 4,0% | 73% | 9,4% | 87% |
| **Média de 6 modelos** | **3,3%** | **90%** | **7,9%** | **96%** |

O «melhor modelo» era escolhido com poucos anos, e o vencedor mudava muitas vezes por acaso: com mais um ano de dados, mudou em 42% dos casos. Os modelos que mais ganhavam eram também os que mais reagiam a pequenas mudanças nos dados, e os que tinham intervalos mais estreitos.

Dois modelos ficaram de fora da média: o passeio aleatório e o modelo ingénuo. Ambos tratam a variação de cada ano como permanente, o que nos municípios pequenos tornava os intervalos tão largos que deixavam de ser úteis.

**O que significa para si.**

- Os intervalos são cerca de duas vezes mais largos do que eram. É a incerteza real a alguns anos de distância.
- A projecção muda menos quando chegam dados novos.
- A escolha do «melhor modelo» continua disponível nas opções avançadas.

### O valor central é o valor esperado

**O que fazemos.** Os modelos trabalham sobre o logaritmo da taxa, e o resultado é convertido de volta com uma correcção, para que o valor central seja o valor esperado e não a mediana.

**Porquê.** O logaritmo impede taxas negativas e trata as variações como proporcionais. Sem a correcção, a projecção ficaria sistematicamente abaixo do esperado, cerca de 5% a 20 anos.

### Uma projecção não é uma meta

**O que fazemos.** A aplicação mostra a projecção com um **intervalo de previsão**.

**Porquê.** O intervalo de previsão diz entre que valores é plausível que um ano futuro venha a ficar. Não é um intervalo de confiança, que mede a incerteza de um valor já observado, e é mais largo do que ele. A projecção diz para onde a série vai se o passado se repetir.

**O que significa para si.** Uma meta é uma decisão: pode usar a projecção como referência, mas não a substitui.

---

## 7. Como se verifica que os números estão certos

**O que fazemos.** Há quatro tipos de verificação:

- **Testes automáticos.** Correm a cada alteração do código.
- **Verificação de regressão.** Corre depois de cada actualização dos dados: recalcula 119 valores validados contra o INE e avisa se algum mudar.
- **Verificação de invariantes.** Também corre depois de cada actualização: calcula os indicadores em muitas áreas, anos e opções, e procura o que nunca deveria acontecer, como um valor em falta onde há dados, Portugal que não bate com a soma dos municípios, ou uma série que salta sem razão.
- **Auditoria das fontes** (`tools/source_audit.R`). Lê os dados tal como o INE os publica e compara cada conjunto com outro que meça o mesmo: municípios com Portugal, sexos com o total, idades com o total, freguesias com municípios, sub-regiões com os seus municípios, semanas com o ano. Foi assim que se encontrou, por exemplo, o erro de 2011 no Alto Alentejo.

**Porquê.** Os erros mais perigosos são os que produzem números plausíveis. Estas verificações já encontraram vários, tanto nossos como nas fontes, e cada um ficou documentado.

**O que significa para si.** Se encontrar um valor estranho, vale a pena comunicá-lo: pode ser um caso que as verificações ainda não cobrem.
