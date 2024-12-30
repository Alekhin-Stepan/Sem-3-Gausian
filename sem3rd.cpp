//112-Alekhin-Stepan-Gaussians

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <ctime>
#include <chrono>
#include <random>          // Для std::mt19937 и std::random_device
#include <iomanip>



//Классы для логирования
class Config {
public:
    int fieldWidth;
    int fieldHeight;
    double defaultX, defaultY, defaultSx, defaultSy, defaultH;
    std::string logFileNameInterface;
    std::string logFileNameControl;
    bool loggingInterfaceEnabled;
    bool loggingControlEnabled;

    Config(const std::string& filename) {
        
        std::ifstream configFile(filename);
        if (!configFile.is_open()) {
            std::cerr << "Failed to open config file." << std::endl;
            return;
        }

     std::string key;
while (configFile >> key) { // Считываем ключи
    if (key == "fieldWidth") configFile >> fieldWidth;
    else if (key == "fieldHeight") configFile >> fieldHeight;
    else if (key == "defaultX") configFile >> defaultX;
    else if (key == "defaultY") configFile >> defaultY;
    else if (key == "defaultSx") configFile >> defaultSx;
    else if (key == "defaultSy") configFile >> defaultSy;
    else if (key == "defaultH") configFile >> defaultH;
    else if (key == "logFileNameInterface") configFile >> logFileNameInterface;
    else if (key == "logFileNameControl") configFile >> logFileNameControl;
    else if (key == "loggingInterfaceEnabled") configFile >> std::boolalpha >> loggingInterfaceEnabled;
    else if (key == "loggingControlEnabled") configFile >> std::boolalpha >> loggingControlEnabled;
}

        configFile.close();
    }
};

class Logger {
private:
    std::ofstream logFile;

public:
    Logger(const std::string& fileName) {
        if (!logFile.is_open()) {
            logFile.open(fileName, std::ios::out | std::ios::app);
        }
        
    }

    ~Logger() {
        if (logFile.is_open()) {
            logFile.close();
        }
    }

    void logMessage(const std::string& message, bool b) {
        if (logFile.is_open() && b) {
            auto now = std::chrono::system_clock::now();
            std::time_t now_c = std::chrono::system_clock::to_time_t(now);
            std::stringstream timeStamp;
            timeStamp << std::put_time(std::localtime(&now_c), "%Y-%m-%d %H:%M:%S");
            logFile << "[" << timeStamp.str() << "] " << message << std::endl;
        }
    }
    
};
    
//Классы для построения гаусса
class Gaus {
public:
    double h, x0, y0, sx, sy;

    Gaus(double h, double x0, double y0, double sx, double sy)
        : h(h), x0(x0), y0(y0), sx(sx), sy(sy) {}
};

class Pole {
public:
    std::vector<std::vector<double>> field;

    Pole(int A, int B) {
        field.resize(A, std::vector<double>(B, 0));
    }
};

class Component {
public:
    std::vector<std::vector<double>> componenta;
    Component(const std::vector<std::vector<double>>& inputComponenta) : componenta(inputComponenta) {}
    
    Component(int A, int B) {
        componenta.resize(A, std::vector<double>(B, 0));
    }
    
    // Метод для получения значений. Например, возвращаем первый вектор
    std::vector<double> getValues() const {
        if (!componenta.empty()) {
            return componenta[0]; // Возвращаем первый вектор
        }
        return {}; // Возвращаем пустой вектор, если componenta пуст
    }

    // Дополнительный метод для получения всех значений в одномерном векторе, если это необходимо
    std::vector<double> getAllValues() const {
        std::vector<double> allValues;
        for (const auto& vec : componenta) {
            allValues.insert(allValues.end(), vec.begin(), vec.end());
        }
        return allValues;
    }
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class GaussBuilder {
   public:
   
       void addgauss(double h, double x0, double y0, double sigma_x, double sigma_y, std::vector<Gaus>& gaussi) {
        gaussi.emplace_back(h, x0, y0, sigma_x, sigma_y);
        //logger.logMessage("Added gauss", b);
    }
    
    void init(int A, int B, Pole*& p) { 
        p = new Pole(A, B);
        //logger.logMessage("Added field", b);
    }
    
    void generate(Pole*& p, std::vector<Gaus>& gaussi) {
    if (p == nullptr) {
       std::cout << "Pole not initialized!" << std::endl;
       return;
   }
   
        double value; 
        for (const auto& g : gaussi) {
            for (long unsigned int x = 0; x < p->field[0].size(); ++x) {
                for (long unsigned int y = 0; y < p->field.size(); ++y) {
                    value = g.h * exp(-((pow((x - g.x0) / g.sx, 2)) + (pow((y - g.y0) / g.sy, 2))) / 2); 
                    p->field[y][x] += value; 
                    if (p->field[y][x] > 255) {
                        p->field[y][x] = 255;
                    }
                    if (p->field[y][x] < 0) {
                        p->field[y][x] = 0;
                    }
                }
            }
        }
    }
};
    
class BmpHandler {
   public:

       void bmp_write(const std::vector<std::vector<double>>& pixelMatrix, const std::string& filename) {
    int width = pixelMatrix[0].size();
    int height = pixelMatrix.size();
    int padding = (4 - (width * 3) % 4) % 4; // Padding for alignment to 4 bytes
    std::ofstream bmpFile(filename, std::ios::binary);
    if (!bmpFile) {
        std::cerr << "Failed to create BMP file." << std::endl;
        //logger.logMessage("Failed to create BMP file.", b);
        return;
    }
    // Write BMP header
    unsigned char bmpHeader[54] = {
        'B', 'M', // Identifier
        0, 0, 0, 0, // Size of file (will be set later)
        0, 0, 0, 0, // Reserved
        54, 0, 0, 0, // Header size
        40, 0, 0, 0, // Info header size
        0, 0, 0, 0, // Width (will be set later)
        0, 0, 0, 0, // Height (will be set later)
        1, 0, // Number of color planes
        24, 0, // Bits per pixel
        0, 0, 0, 0, // Compression
        0, 0, 0, 0, // Image size (will be set later)
        0x13, 0x0B, 0, 0, // Horizontal resolution
        0x13, 0x0B, 0, 0, // Vertical resolution
        0, 0, 0, 0, // Number of colors in palette
        0, 0, 0, 0  // Important colors
    };
    // Set width and height in header
    bmpHeader[18] = (width & 0xFF);
    bmpHeader[19] = (width >> 8) & 0xFF;
    bmpHeader[20] = (width >> 16) & 0xFF;
    bmpHeader[21] = (width >> 24) & 0xFF;
    bmpHeader[22] = (height & 0xFF);
    bmpHeader[23] = (height >> 8) & 0xFF;
    bmpHeader[24] = (height >> 16) & 0xFF;
    bmpHeader[25] = (height >> 24) & 0xFF;
    // Write header
    bmpFile.write(reinterpret_cast<char*>(bmpHeader), 54);
    // Write pixel data
    for (int y = height - 1; y >= 0; --y) { // BMP stores pixels bottom-to-top
        for (int x = 0; x < width; ++x) {
            unsigned char color = 255 - static_cast<unsigned char>(pixelMatrix[y][x]); // Color
            bmpFile.put(color); // B
            bmpFile.put(color); // G
            bmpFile.put(color); // R
        }
        // Add padding
        for (int p = 0; p < padding; ++p) {
            bmpFile.put(0);
        }
    }
    bmpFile.close();
}
    
   void bmp_read(GaussBuilder& gaussBuilder, const std::string &filename, std::vector<std::vector<double>> &pixelMatrix, Pole*& p) {
    std::ifstream bmpFile(filename, std::ios::binary);
    if (!bmpFile) {
        std::cerr << "Failed to open BMP file." << std::endl;
        //logger.logMessage("Failed to open BMP file.", b);
        return;
    }
    // Читаем заголовок BMP
    unsigned char header[54];
    bmpFile.read(reinterpret_cast<char*>(header), 54);
    // Получаем ширину и высоту изображения
    int width = header[18] | (header[19] << 8) | (header[20] << 16) | (header[21] << 24);
    int height = header[22] | (header[23] << 8) | (header[24] << 16) | (header[25] << 24);
    // Инициализируем новое поле
    gaussBuilder.init(height, width, p); // Заметь, что BMP хранит данные от нижней строки к верхней.
    
    // Читаем данные пикселей
    for (int y = height - 1; y >= 0; --y) { // BMP хранит данные снизу вверх
        for (int x = 0; x < width; ++x) {
            unsigned char color = bmpFile.get(); // Читаем B
            bmpFile.get(); // Читаем G
            bmpFile.get(); // Читаем R
            double value = 255 - color; // Цвет в высоту
            pixelMatrix[y][x] = value; // Обновляем матрицу значений
            //Значения в gaussi некорректные
        }
        bmpFile.ignore((4 - (width * 3) % 4) % 4); // Пропускаем паддинг
    }
    bmpFile.close();
  }
};
   
class GnuplotInterface {
   public:
   
       void gnuplot(Pole*& p) {
       if (p == nullptr) {
       std::cout << "Pole not initialized!" << std::endl;
       return;
   }
        int rows = p->field.size(); 
        int cols = p->field[0].size(); 
        // Open a pipe to gnuplot 
        FILE* gnuplotPipe = popen("gnuplot -p", "w"); 
        if (!gnuplotPipe) { 
            std::cerr << "Could not open pipe to gnuplot." << std::endl; 
            //logger.logMessage("Could not open pipe to gnuplot.", b);
            return; 
        }
        // Send gnuplot commands 
        fprintf(gnuplotPipe, "set contour base\n"); 
        fprintf(gnuplotPipe, "set view 60,30\n"); 
        fprintf(gnuplotPipe, "set xrange [0:%d]\n", cols - 1); 
        fprintf(gnuplotPipe, "set yrange [0:%d]\n", rows - 1); 
        fprintf(gnuplotPipe, "set zrange [0:255]\n"); // Set z range based on 0 and 255
        fprintf(gnuplotPipe, "set terminal png\n"); 
        fprintf(gnuplotPipe, "set output 'landscape.png'\n"); 
        fprintf(gnuplotPipe, "splot '-' with lines\n"); 
        // Write data directly to gnuplot 
        for (int y = 0; y < rows; ++y) { 
            for (int x = 0; x < cols; ++x) { 
                fprintf(gnuplotPipe, "%d %d %f\n", x, y, p->field[y][x]); 
            } 
            fprintf(gnuplotPipe, "\n"); // Newline to separate rows 
        } 
        fprintf(gnuplotPipe, "EOF\n"); // Close the pipe 
        pclose(gnuplotPipe);
    }
};
   
class ComponentCalculator {
   public:
    int count = 0;//Для шума
    int Noise = congig.Noise; //шум из конфига
       
    int incrementAndCollect(std::vector<std::vector<double>>& componenta, std::vector<std::vector<double>> &CopyPole, int x, int y, int i) {

        if (x < 1 || y < 1 || x > (int) componenta[0].size() - 2 || y > (int) componenta.size() - 2 || CopyPole[y][x] < 250) return -1;

        if (CopyPole[y][x] >= 255 && CopyPole[y][x] <= 255) {
            CopyPole[y][x] = 0; // Пометить как посещенное
            count = count < i + 1 ? i + 1 : count;
            componenta[y][x] = 255; // Увеличить значение в Componenta
            incrementAndCollect(componenta, CopyPole, x + 1, y, i + 1);
            incrementAndCollect(componenta, CopyPole, x - 1, y, i + 1);
            incrementAndCollect(componenta, CopyPole, x, y + 1, i + 1);
            incrementAndCollect(componenta, CopyPole, x, y - 1, i + 1);
        }
        
    return count;
}
       void bin(BmpHandler& bmpHandler, std::vector<std::vector<double>> &CopyPole, int slise, Pole*& p) {
       if (p == nullptr) {
       std::cerr << "Pole not initialized!" << std::endl;
       return;
   }
       CopyPole = p->field;//Копия
        
            for (int x = 0; x < (int) p->field[0].size(); ++x) {
                for (int y = 0; y < (int) p->field.size(); ++y) {
                   CopyPole[y][x] = p->field[y][x] > slise ? 255 : 0;
                } 
            }
        
        bmpHandler.bmp_write(CopyPole, "slise.bmp");
        //loggerinterface.logMessage("Created BMP file.", b);
    }
    
    void wave(std::vector<Component>& componenti, std::vector<std::vector<double>> &CopyPole, Pole*& p) {
    if (p == nullptr) {
       std::cerr << "Pole not initialized!" << std::endl;
       return;
   }

            for (int y = 0; y < (int) p->field.size(); ++y) {
                for (int x = 0; x < (int) p->field[y].size(); ++x) {
                    if (CopyPole[y][x] <= 255 && CopyPole[y][x] >= 255) {
                        count = 0;
                        Component Componenta(p->field.size(), p->field[0].size());
                        if (incrementAndCollect(Componenta.componenta, CopyPole, x, y, 1) >= Noise) {
                         componenti.emplace_back(Componenta);
                       }
                   }
               }
           }

        loggerinterface.logMessage("Wave used, amount component = " + std::to_string(componenti.size()), b);
   }
};




class Control {
private:
    bool b = true; 
public:
    //Для логирования
    Config& config;
    Logger& loggercontrol;
    //Для функций
    std::vector<std::vector<double>> CopyPole;
    std::vector<Gaus> gaussi;
    std::vector<Component> componenti;
    GaussBuilder gaussBuilder;
    BmpHandler bmpHandler;
    GnuplotInterface gnuplotInterface;
    ComponentCalculator componentCalculator;
    Pole* p = nullptr; // Pointer to Pole
    
     
    Control(Config& cfg, Logger& log) : config(cfg), loggercontrol(log) {
        if (config.loggingControlEnabled) {
            
                loggercontrol.logMessage("Logging Control is enabled.", b);
                std::cout << "Logging Control is enabled." << std::endl;
            } else {
            loggercontrol.logMessage("Logging Control is disabled.", b);
            std::cout << "Logging Control is disabled." << std::endl;
            b = false;
        }      
    }
    
    ~Control() {
        delete p;
    }
    
void Dispetcher(const std::string& command, int A, int B, double h, double x0, double y0, double sigma_x, double sigma_y, 
std::vector<std::vector<double>>& pixelMatrix, const std::string& filename, int slise, int klaster) {
   
    if (command == "init") {
    gaussBuilder.init(A, B, p);
    }

    if (command == "g") {
    gaussBuilder.addgauss(h, x0, y0, sigma_x, sigma_y, gaussi);
    }
    
    if (command == "generate") {
    gaussBuilder.generate(p, gaussi);
    }
    
    if (command == "gnuplot") {
    gnuplotInterface.gnuplot(p);
    }

    if (command == "bmp_write") {
    bmpHandler.bmp_write(pixelMatrix, filename);
    }

    if (command == "bmp_read") {
    bmpHandler.bmp_read(gaussBuilder, filename, pixelMatrix, p);
    }

    if (command == "bin") {
    componentCalculator.bin(bmpHandler, CopyPole, slise, p);
    componentCalculator.wave(componenti, CopyPole, p);
    }
    
    
        }
};
   //Деструктор!!!!!!!!!!!


class Interface {
private:
bool b = true;

public:
    Control& c;
    Config& config;
    Logger& loggerinterface;
    
    Interface(Config& cfg, Logger& log, Control& c) : config(cfg), loggerinterface(log), c(c){       
        if (config.loggingInterfaceEnabled) {
            
                loggerinterface.logMessage("Logging Interface is enabled.", b);
                std::cout << "Logging Interface is enabled." << std::endl;
        } else {
        loggerinterface.logMessage("Logging Interface is disabled.", b);
        std::cout << "Logging Interface is disabled." << std::endl;
        b = false;
      }        
  }
  
  // Получаем размеры поля из конфигурации
    int A = config.fieldWidth;  // Размеры из конфигурации
    int B = config.fieldHeight;
    int slice = 0;
    int k = 0;//cluster
    
    void print() {
        double x = 0, y = 0, sx = 0, sy = 0, h = 0;
        std::string s;
        bool a;
        std::string filename;
        std::ifstream file;
        std::cout << "112-Alekhin-Gaussians-3rd-sem.\nHello, dear user, this program builds Gaussians.\nEnter commands from a text file (PRESS 0) or from the keyboard (PRESS 1)?" << std::endl;
        std::cin >> a;
        loggerinterface.logMessage("User chose input method: " + std::to_string(a), b);
        if (a == 0) {
            std::cout << "You will enter commands from a text file.\nEnter filename:" << std::endl;
            std::cin >> filename;
            loggerinterface.logMessage("Reading commands from file: " + filename, b);
            file.open(filename);
            if (!file) {
                std::cout << "File not found" << std::endl;
                loggerinterface.logMessage("Error: File not found.", b);
                return;
            }
        } else {
            std::cout << "You will enter commands from the keyboard" << std::endl;
            loggerinterface.logMessage("User chose to input commands from the keyboard.", b);
        }
        if (a == 0) {
            int n = 0;
            while (file >> s) {
                loggerinterface.logMessage("Received command: " + s, b);
         if (s == "help") {
            // Открываем файл help.txt для записи
            std::ofstream helpFile("help.txt");
            if (helpFile.is_open()) {
                helpFile << "Available commands:\n";
                helpFile << "help - создание файла с пояснением команд\n";
                helpFile << "init - инициализация поля. Первая команда!!!\n";
                helpFile << "g(x,y,sx,sy,h) - создает гаусс с 5 параметрами\n";
                helpFile << "generate - складывает гауссы\n";
                helpFile << "gnuplot - рисует картинку в gnuplot\n";
                helpFile << "bmp_write - создает черно-белую серую картинку bmp\n";
                helpFile << "bmp_read - чтение bmp файла и инициализация поля новыми размерами\n";
                helpFile << "bin - срез: все что выше черное, все что ниже белое, считается число компонент, каждая компонента запомнена\n";
                helpFile << "k_means - выделение k кластеров\n";
                helpFile << "end - это конец программы\n";
                helpFile.close(); // Закрываем файл
                loggerinterface.logMessage("Help file created successfully.", b);
            } else {
                loggerinterface.logMessage("Failed to create help file.", b);
            }
        }
   else if (s == "init") {
    // Проверяем, была ли уже вызвана команда init
    if (n != 0) {
        std::cout << "The init command has already been called.\nError\n";
        loggerinterface.logMessage("Error: Multiple init commands.", b);
        return;
    }
     n = 1;
    // Логируем и инициализируем поле
    loggerinterface.logMessage("Initializing field with size: " + std::to_string(A) + " x " + std::to_string(B), b);
    c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, "output.bmp", slice, k);
    loggerinterface.logMessage("Field initialized.", b);
} else if (n != 1) {
                     std::cout << "The init command was not use.\nError\n";
                     loggerinterface.logMessage("Error: The init command was not use.", b);
                     return;
             }
 else if (s == "g") {

                // Читаем параметры из файла
                file >> x >> y >> sx >> sy >> h;

                // Проверяем, были ли параметры указаны
                if (file.fail()) {
                    // Если не удалось прочитать, значит, используем значения по умолчанию
                    if (file.eof()) { // Желательно проверить конец файла, чтобы избежать повторного чтения
                        x = config.defaultX;
                        y = config.defaultY;
                        sx = config.defaultSx;
                        sy = config.defaultSy;
                        h = config.defaultH;
                    } else {
                        // Если произошла ошибка чтения или недостаточно параметров
                        // Нужно сбросить флаг ошибки для дальнейшего чтения
                        file.clear();

                        // Проверяем, какие параметры были прочитаны
                        if (!(file >> x)) x = config.defaultX; // Если не удалось получить x, берем по умолчанию
                        if (!(file >> y)) y = config.defaultY; // Если не удалось получить y, берем по умолчанию
                        if (!(file >> sx)) sx = config.defaultSx; // ...
                        if (!(file >> sy)) sy = config.defaultSy; // ...
                        if (!(file >> h)) h = config.defaultH; // ...
                    }
                }

                loggerinterface.logMessage("Adding Gaussian: x=" + std::to_string(x) +
                    ", y=" + std::to_string(y) + 
                    ", sx=" + std::to_string(sx) + 
                    ", sy=" + std::to_string(sy) + 
                    ", h=" + std::to_string(h), b);
                c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, "output.bmp", slice, k);
            } else if (s == "generate") {
                    c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, "output.bmp", slice, k);
                    loggerinterface.logMessage("Generated values in the field.", b);
                } else if (s == "gnuplot") {
                    c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, "output.bmp", slice, k);
                    loggerinterface.logMessage("Called gnuplot.", b);
                } else if (s == "bmp_write") {
                    c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, "output.bmp", slice, k);
                    loggerinterface.logMessage("Created BMP file.", b);
                } else if (s == "bmp_read") {
                    file >> filename; // Чтение имени файла для bmp_read
                    c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, filename, slice, k);
                    loggerinterface.logMessage("Read BMP file: " + filename, b);
                } else if (s == "bin") {
                    file >> slice;
                    c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, "output.bmp", slice, k);
                    loggerinterface.logMessage("Slice applied: slice=" + std::to_string(slice), b);
                    loggerinterface.logMessage("Wave will be used", b);
                }
            }
        } else {
            int n = 0;
            while (true) {
                std::cout << "Enter command (help, init, g, generate, gnuplot, bmp_write, bmp_read, bin, end):";
                std::cin >> s;
                std::cout << "\n";
                loggerinterface.logMessage("Received command: " + s, b);
        if (s == "help") {
            // Открываем файл help.txt для записи
            std::ofstream helpFile("help.txt");
            if (helpFile.is_open()) {
                helpFile << "Available commands:\n";
                helpFile << "help - создание файла с пояснением команд\n";
                helpFile << "init - инициализация поля. Первая команда!!!\n";
                helpFile << "g(x,y,sx,sy,h) - создает гаусс с 5 параметрами\n";
                helpFile << "generate - складывает гауссы\n";
                helpFile << "gnuplot - рисует картинку в gnuplot\n";
                helpFile << "bmp_write - создает черно-белую серую картинку bmp\n";
                helpFile << "bmp_read - чтение bmp файла и инициализация поля новыми размерами\n";
                helpFile << "bin - срез: все что выше черное, все что ниже белое, считается число компонент, каждая компонента запомнена\n";
                helpFile << "k_means - выделение k кластеров\n";
                helpFile << "end - это конец программы\n";
                helpFile.close(); // Закрываем файл
                loggerinterface.logMessage("Help file created successfully.", b);
            } else {
                loggerinterface.logMessage("Failed to create help file.", b);
            }
        }
    if (s == "init") {
    // Проверяем, была ли уже вызвана команда init
    if (n != 0) {
        std::cout << "The init command has already started.\nError\n";
        loggerinterface.logMessage("Error: Multiple init commands.", b);
        return;
    }
    n = 1; // Устанавливаем флаг, что команда init была вызвана
    
    // Логируем и инициализируем поле
    loggerinterface.logMessage("Initializing field with size: " + std::to_string(A) + " x " + std::to_string(B), b);
    c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, "output.bmp", slice, k);
    loggerinterface.logMessage("Field initialized.", b);
}
  if (n != 1) {
                     std::cout << "The init command was not use.\nError\n";
                     loggerinterface.logMessage("Error: The init command was not use.", b);
                     return;
             }
             
  if (s == "g") {
    std::string input;
    std::getline(std::cin, input);

    std::istringstream inputStream(input);

    // Инициализируем переменные значениями по умолчанию
    x = config.defaultX;
    y = config.defaultY;
    sx = config.defaultSx;
    sy = config.defaultSy;
    h = config.defaultH;

    // Читаем введенные данные
    if (!(inputStream >> x)) {
        std::cout << "The default value for x is used: " << config.defaultX << std::endl;
    }
    if (!(inputStream >> y)) {
        std::cout << "The default value for y is used: " << config.defaultY << std::endl;
    }
    if (!(inputStream >> sx)) {
        std::cout << "The default value for sx is used: " << config.defaultSx << std::endl;
    }
    if (!(inputStream >> sy)) {
        std::cout << "The default value for sy is used: " << config.defaultSy << std::endl;
    }
    if (!(inputStream >> h)) {
        std::cout << "The default value for h is used: " << config.defaultH << std::endl;
    }

                loggerinterface.logMessage("Adding Gaussian: x=" + std::to_string(x) +
                    ", y=" + std::to_string(y) + 
                    ", sx=" + std::to_string(sx) + 
                    ", sy=" + std::to_string(sy) + 
                    ", h=" + std::to_string(h), b);
               c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, "output.bmp", slice, k);
            }

                  if (s == "generate") {
                    c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, "output.bmp", slice, k);
                    std::cout << "Generated values in the field." << std::endl;
                    loggerinterface.logMessage("Generated values in the field.", b);
                } if (s == "gnuplot") {
                    c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, "output.bmp", slice, k);
                    std::cout << "Called gnuplot." << std::endl;
                    loggerinterface.logMessage("Called gnuplot.", b);
                } if (s == "bmp_write") {
                    c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, "output.bmp", slice, k);
                    std::cout << "Created BMP file." << std::endl;
                    loggerinterface.logMessage("Created BMP file.", b);
                } if (s == "bmp_read") {
                    std::cout << "Enter the filename to read:" << std::endl;
                    std::cin >> filename;
                    c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, filename, slice, k);
                    std::cout << "Read BMP file: " + filename << std::endl;
                    loggerinterface.logMessage("Read BMP file: " + filename, b);
                }  if (s == "end") {
                    std::cout << "Ending the program" << std::endl;
                    loggerinterface.logMessage("Ending the program.", b);
                    break;
                }
                   if (s == "bin") {
                     std::cout << "Enter slice level:" << std::endl;
                     std::cin >> slice;
                     c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, "output.bmp", slice, k);
                     loggerinterface.logMessage("Slice applied: slice=" + std::to_string(slice), b);
                     std::cout << "Slice applied: slice=" << slice << std::endl;
                     loggerinterface.logMessage("Wave will be used", b);
                }
                    
            }

        if (file.is_open()) {
            file.close();
            loggerinterface.logMessage("Closed input file.", b);
        }
     }
   }  
};

int main() {
    //Логирование
    Config config("config.txt");
    Logger loggerinterface(config.logFileNameInterface);
    Logger loggercontrol(config.logFileNameControl);
    // Создаем интерфейс
    Control c(config, loggercontrol);
    Interface i(config, loggerinterface, c);
    
    // Вызываем метод print() интерфейса
    i.print();

    return 0;
}
