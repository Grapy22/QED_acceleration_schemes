FC = gfortran
FFLAGS = -Wall -O3 -march=native -funroll-loops -ftree-vectorize
CUBA_DIR = /home/dmitry/Cuba-4.2.2
LIBS = -L$(CUBA_DIR) -lcubaq -lm

# Явно перечисляем ВСЕ исходные файлы
SOURCES = $(wildcard *.f90 *.f)

# Создаем список объектных файлов
OBJECTS = $(SOURCES:.f90=.o)
OBJECTS := $(OBJECTS:.f=.o)

TARGET = main

# Цель по умолчанию зависит от объектных файлов
all: $(TARGET)

# Сборка программы из объектных файлов
$(TARGET): $(OBJECTS)
	$(FC) $(OBJECTS) -o $(TARGET) $(LIBS)
	@echo "✅ Сборка завершена!"

# Правило для .f90 -> .o
%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@
	@echo "Компилирую $<"

# Правило для .f -> .o
%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@
	@echo "Компилирую $<"

# Запуск
run: $(TARGET)
	./$(TARGET)

# Очистка
clean:
	rm -f $(TARGET) *.o *.mod
	@echo "🧹 Очистка завершена!"

# Полная пересборка
rebuild: clean all
