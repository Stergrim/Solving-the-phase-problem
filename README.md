# Решение фазовой проблемы <br> (Solving the phase problem)

В данном репозитории представлена программная реализация решения фазовой проблемы путём многомерной оптимизации волнового фронта по картине распределения интенсивности.

Подход заключается в том, чтобы путём сравнения изображений восстановить волновой фронт излучения зная только распределение интенсивности. Для минимизации многомерной функции используется алгоритм ***БФГШ***.

## Каталоги в этом репозитории

>**demos**: папка, содержащая демонстрацию работы предложенного решения <br>
>**matlab**: папка, содержащая код программы

## Дипломная работа

Бакалаврская работа:<br>

Находится в папке `demos` под названием **DiplomaBachelor.pdf** <br>
**Google Disk:** https://drive.google.com/file/d/10HKgejvoL49USyzA6oxwIyHwyJ5L9aMx/view?usp=sharing <br>
**ЯндексДиск:** https://disk.yandex.ru/i/R19TsmHZzovZfg

## Тестирование

Подробное описание запуска программы и пример изображения представлен в папке `demos` в файле [README](https://github.com/Stergrim/Solving-the-phase-problem/blob/main/demos/README.md)

## Результаты проверки решения на синтезированных изображениях

Пример исходного изображения:

<figure>
<img src="https://github.com/Stergrim/Solving-the-phase-problem/blob/main/demos/GroundTruthModelVisual.png" width="250"/>
</figure>

Изображение с которого начинается расчёт и восстановленное изображение:

<p float="left">
<img src="https://github.com/Stergrim/Solving-the-phase-problem/blob/main/demos/StartModelVisual.png" width="250" />
<img src="https://github.com/Stergrim/Solving-the-phase-problem/blob/main/demos/RestoredModelVisual.png" width="250" /> 
</p>

Матраца разности изображений, максимум на изображении равен **4**, и восстановленный волновой фронт:

<figure>
<img src="https://github.com/Stergrim/Solving-the-phase-problem/blob/main/demos/DifferenceModelVisual.png" width="250"/>     <img src="https://github.com/Stergrim/Solving-the-phase-problem/blob/main/demos/WaveFrontModelVisual.png" width="250"/>
</figure>

СКО изображений **1,18** пикселя, погрешность восстановления волнового фронта **0,12%**.

## Результаты решения на реальных изображениях

Стенд для получения изображений:<br>
1 – источник излучения; 2 – конденсор; 3 – тест-объект; 4 – объектив коллиматора; 5 – исследуемый объектив; 6 – система регистрации изображения; 7 – поворотный узел; 8 – светозащитный чехол.

<figure>
<img src="https://github.com/Stergrim/Solving-the-phase-problem/blob/main/demos/ExperimentalSetup.png" width="700"/>
</figure>

Полученные со стенда изображения, под углом в **5** и **10** градусов к оптической оси исследуемого объектива:

<figure>
<img src="https://github.com/Stergrim/Solving-the-phase-problem/blob/main/demos/ExperimentReal.png" width="550"/>
</figure>

Теоретические изображения, рассчитанные в ***Zemax***:

<figure>
<img src="https://github.com/Stergrim/Solving-the-phase-problem/blob/main/demos/ExperimentTheor.png" width="550"/>
</figure>

Теоретические изображения использовались в качестве начального приближения и в результате расчёта **СКО** составил **1,09** и **1,17** пикселя. <br>

Восстановленные волновые фронты:

<figure>
<img src="https://github.com/Stergrim/Solving-the-phase-problem/blob/main/demos/ExperimentWaveFronts.png" width="700"/>
</figure>

## Замечания

К решению этой задачи я подходил, совершенно ничего не зная о программировании, и поэтому весь код написан довольно неуклюже и просто, без использования подходов *ООП* и *ФП*.
Минусами моей реализации являются:
1.	Много параметров, которые необходимо задать для начала расчёта и из-за отсутствия наглядности (визуализации) можно запутаться в последовательности действий для запуска.
2.	Время расчёта. Как минимум можно значительно сократить время расчёта, если в функции `DirectTask` рассчитывать распределение интенсивности не напрямую, а через преобразование Фурье. <br>

Буду рад аргументированной критике моего решения и советам по улучшению программы.
