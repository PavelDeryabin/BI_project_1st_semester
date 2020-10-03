# Осенний проект в Институте Биоинформатики

## Связь

Почта: deryabin.pav@gmail.com 
Телеграм: @deryabin_pavel
Телефон: +7-911-712-40-88

# Структура репозитория /

## Полезные работы по теме /articles

- Клеточное старение (senescence) и таргетная элиминация стареющих клеток (senolysis) **/articles/basic_background**
- Пилотные работы по использованию сердечных гликозидов (cardiac glycosides) как сенолитиков **/articles/cg_senolysis**

## Данные /data

- Сырые данные **/data/raw_data**
    - Наши данные - эндометриальные МСК (eMSC), 4 контрольных образца и 4 образца состаренных клеток (преждевременное стресс-индуцированное старение, Н2О2, 1 час, 200 мкМ, 2 недели после стресса). *На этих клетках в наших экспериментах мы наблюдаем отсутствие сенолитического эффекта при использовании сердечного гликозида оуабаина.* **/data/raw_data/emsc**
    - [GSE122081](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122081) - IMR90 ER:RAS, 3 контрольных образца и 3 образца состаренных клеток (преждевременное онкоген-индуцированное старение, индуцибельная оверэкспрессия RAS, неделя после индукции). *Это данные непосредственно тех клеток, на которых сенолизиз оуабаином имел место, согласно пилотной работе [1]. Ввиду отсутствия такой линии в лаборатории, экспериментально эффект на этих клетках мы не валидировали.* **/data/raw_data/imr90**
    - [GSE102639](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102639) - A549, 2 контрольных образца и 2 образца состаренных клеток (преждевременное стресс-индуцированное старение, этопозид, 2 мкМ, неделя после стресса). *Это данные статьи [2], где не оценивалось сенолитическое действие оуабаина, но использовалась модель старения A549 с практически такими же условиями индукции старения, какие использовались в пилотных работах [1, 3]. При воспроизведении экспериментов пилотных работ в наших условиях, сенолитический эффект оуабаина на этих клетках воспроизводится.* **/data/raw_data/a549**
    1. [Guerrero A. et al. Cardiac glycosides are broad-spectrum senolytics // Nature metabolism. – 2019. – Т. 1. – №. 11. – С. 1074-1088.](https://www.nature.com/articles/s42255-019-0122-z)
    2. [Wang L. et al. High-throughput functional genetic and compound screens identify targets for senescence induction in cancer // Cell reports. – 2017. – Т. 21. – №. 3. – С. 773-783.](https://www.cell.com/cell-reports/fulltext/S2211-1247(17)31390-6?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2211124717313906%3Fshowall%3Dtrue)
    3. [Triana-Martínez F. et al. Identification and characterization of Cardiac Glycosides as senolytic compounds // Nature communications. – 2019. – Т. 10. – №. 1. – С. 1-12.](https://www.nature.com/articles/s41467-019-12888-x)

    *в основной директории **/data/raw_data/README** с описанием загрузки данных и соответствующие bash-скрипты

    **статьи 1,2,3 лежат в репозитории  **/articles/cg_senolysis**

- Контроль качества ридов **/data/qc**
    - fastQC отчеты по сырым данным **/data/qc/fastqc_of_raw_data**
    - fastQC отчеты по обработанным данным **/data/qc/fastqc_of_processed_data**
    - **/data/qc/README** с описанием примененных операций фильтрации и тримминга, и соответствующие bash-скрипты

- Выравнивание/квантификация **/data/salmon**
    - Индекс референсного транскриптома **/data/salmon/index**
    - Результаты меппинга **/data/salmon**
    - **/data/salmon/README** с описанием SA квантификации и соответствующими bash-скриптами

- tximeta квантификация, EDA и DE в DESeq2, биологическая интерпретация **/data/ouab**
    - Дизайн эксперимента **/data/ouab/coldata.csv**
    - Основной R-скрипт работы с данными **/data/ouab/ouab.R** и рабочее окружение **/data/ouab/ouab.RData**
    - Ключевые результаты анализа **/data/ouab/results**