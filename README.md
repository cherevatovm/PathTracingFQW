# Трассировка пути
Проект, реализующий глобальное освещение методом трассировки пути.
Используется структура ускорения поиска пересечений BVH.
Реализована поддержка следующих видов поверхностей:
1. Гладких — диффузной, зеркальной, прозрачной.
2. Шероховатых (на основе микрограней) — отражающей, прозрачной.
Приложение позволяет отрендерить сцену Корнуэльской комнаты с двумя объектами, указав разрешение изображения и количество семплов на пиксель.
3D-модель, а также материал для каждого объекта можно выбрать из заготовленного списка.
После завершения рендеринга изображение может быть сохранено в желаемую директорию с помощью графического интерфейса.


# Path tracing
A project implementing global illumination using the path tracing method.
Utilizes a BVH ray intersection acceleration structure.
Supports the following surface types:
1. Smooth – diffuse, specular, transparent.
2. Rough (microfacet-based) – reflective, transparent.
The application allows to render a Cornell Box scene with two objects, where the user can specify the image resolution and the number of samples per pixel.
A 3D model, as well as a material for each object, can be selected from a list.
After rendering is complete, the image can be saved to a desired directory via the graphical interface.
