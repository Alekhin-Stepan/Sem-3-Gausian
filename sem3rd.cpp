//112-Alekhin-Stepan-Gaussians

#include <iostream>
#include <algorithm>
#include <string>
#include <fstream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <ctime>
#include <chrono>
#include <random>          
#include <iomanip>



//Классы для логирования
class Config {
public:
    int fieldWidth;
    int fieldHeight;
    double defaultX, defaultY, defaultSx, defaultSy, defaultH;
    int Noise;
    double R;
    double pitch, roll;
    std::string logFileNameInterface;
    std::string logFileNameControl;
    bool loggingInterfaceEnabled;
    bool loggingControlEnabled;
    double eps;
    

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
    else if (key == "Noise") configFile >> Noise;
    else if (key == "R") configFile >> R;
    else if (key == "pitch") configFile >> pitch;
    else if (key == "roll") configFile >> roll;
    else if (key == "epsilon") configFile >> eps;
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

struct point {
    double x, y;
    point(double x_ = 0, double y_ = 0) : x(x_), y(y_) {}
    
    bool operator==(const point& other) const {
        return std::fabs(x - other.x) < 0.0001 && std::fabs(y - other.y) < 0.0001;
    }
    
    bool operator!=(const point& other) const {
        return !(*this == other);
    }
};

struct edge {
    point a, b;
    edge(point a_, point b_) : a(a_), b(b_) {}
    
    bool operator==(const edge& other) const {
        return (a == other.a && b == other.b) || (a == other.b && b == other.a);
    }
};

struct triangle {
    point a, b, c;
    triangle(point a_, point b_, point c_) : a(a_), b(b_), c(c_) {}

    // Добавляем оператор сравнения
    bool operator==(const triangle& other) const {
        return (a == other.a && b == other.b && c == other.c) ||
               (a == other.a && b == other.c && c == other.b) ||
               (a == other.b && b == other.a && c == other.c) ||
               (a == other.b && b == other.c && c == other.a) ||
               (a == other.c && b == other.a && c == other.b) ||
               (a == other.c && b == other.b && c == other.a);
    }
    
    point CircumCenter() const  //функция поиска центра описанной окружности
    {
        return point(-0.5*((b.x*b.x-c.x*c.x+b.y*b.y-c.y*c.y)*(a.y-b.y)-(a.x*a.x-b.x*b.x+a.y*a.y-b.y*b.y)*(b.y-c.y))/((a.x-b.x)*(b.y-c.y)-(b.x-c.x)*(a.y-b.y)), 0.5*((b.x*b.x-c.x*c.x+b.y*b.y-c.y*c.y)*(a.x-b.x)-(a.x*a.x-b.x*b.x+a.y*a.y-b.y*b.y)*(b.x-c.x))/((a.x-b.x)*(b.y-c.y)-(b.x-c.x)*(a.y-b.y)));
    }
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class GaussBuilder {
   bool b=true;
   public:
   
       void addgauss(double h, double x0, double y0, double sigma_x, double sigma_y, std::vector<Gaus>& gaussi, Logger &logger, bool b) {
        gaussi.emplace_back(h, x0, y0, sigma_x, sigma_y);
        logger.logMessage("Added gauss", b);
    }
    
    void init(int A, int B, Pole*& p, Logger &logger, bool b) { 
        p = new Pole(A, B);
        logger.logMessage("Added field", b);
    }
    
    void generate(Pole*& p, std::vector<Gaus>& gaussi, Logger &logger, bool b) {
    if (p == nullptr) {
       std::cout << "Pole not initialized!" << std::endl;
       logger.logMessage("Error while generating, field wasn't initialized!", b);
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
            logger.logMessage("Added gaussian into field", b);
        }
    }
};
    
class BmpHandler {
   
   public:

       void bmp_write(const std::vector<std::vector<double>>& pixelMatrix, const std::string& filename, Logger &logger, bool b) {
    int width = pixelMatrix[0].size();
    int height = pixelMatrix.size();
    int padding = (4 - (width * 3) % 4) % 4; // Padding for alignment to 4 bytes
    std::ofstream bmpFile(filename, std::ios::binary);
    if (!bmpFile) {
        std::cerr << "Failed to create BMP file." << std::endl;
        logger.logMessage("Failed to create BMP file.", b);
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
    logger.logMessage("Bmp file created successfully", b);
}
    
   void bmp_read(GaussBuilder& gaussBuilder, const std::string &filename, std::vector<std::vector<double>> &pixelMatrix, Pole*& p, Logger &logger, bool b) {
    std::ifstream bmpFile(filename, std::ios::binary);
    if (!bmpFile) {
        std::cerr << "Failed to open BMP file." << std::endl;
        logger.logMessage("Failed to open BMP file.", b);
        return;
    }
    // Читаем заголовок BMP
    unsigned char header[54];
    bmpFile.read(reinterpret_cast<char*>(header), 54);
    // Получаем ширину и высоту изображения
    int width = header[18] | (header[19] << 8) | (header[20] << 16) | (header[21] << 24);
    int height = header[22] | (header[23] << 8) | (header[24] << 16) | (header[25] << 24);
    // Инициализируем новое поле
    gaussBuilder.init(height, width, p, logger, b); // Заметим, что BMP хранит данные от нижней строки к верхней.
    logger.logMessage("New field initialized.", b);
    
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
   
       void gnuplot(Pole*& p, Logger &logger, bool b) {
       if (p == nullptr) {
       std::cout << "Pole not initialized!" << std::endl;
       logger.logMessage("Error whith gnuplot, field wasn't initialized!", b);
       
       return;
   }
        int rows = p->field.size(); 
        int cols = p->field[0].size(); 
        // Open a pipe to gnuplot 
        FILE* gnuplotPipe = popen("gnuplot -p", "w"); 
        if (!gnuplotPipe) { 
            std::cerr << "Could not open pipe to gnuplot." << std::endl; 
            logger.logMessage("Could not open pipe to gnuplot.", b);
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
   
       
    void incrementAndCollect(std::vector<std::vector<double>>& componenta, std::vector<std::vector<double>> &CopyPole, unsigned short x, unsigned short y, int i) {
        
        

        if (CopyPole[y][x] >= 255 && CopyPole[y][x] <= 255) {
            CopyPole[y][x] = 0; // Пометить как посещенное
            count = count < i + 1 ? i + 1 : count; 
            
            // std::cout << "Added (" << x << ", " << y << ", " << i << ", " << count <<") ";
            
            componenta[y][x] = 255; // Увеличить значение в Componenta
            
            
            if ( x < (int)componenta[0].size()-1 ) { if ( CopyPole[y][x+1] >= 250 ) {
                incrementAndCollect(componenta, CopyPole, x + 1, y, i + 1);
            }}
            if ( x > 1 ) { if ( CopyPole[y][x-1] >= 250 ) {
                incrementAndCollect(componenta, CopyPole, x - 1, y, i + 1);
            }}
            if ( y < (int)componenta.size()-1 ) { if ( CopyPole[y+1][x] >= 250 ) {
                incrementAndCollect(componenta, CopyPole, x, y + 1, i + 1);
            }}
            if ( y > 1 ) { if (CopyPole[y-1][x] >= 250) {
                incrementAndCollect(componenta, CopyPole, x, y - 1, i + 1);
            }}
        }
       
    return;
}
       void bin(BmpHandler& bmpHandler, std::vector<std::vector<double>> &CopyPole, int slise, const std::string &filename, Pole*& p, Logger &logger, bool b) {
       if (p == nullptr) {
       std::cerr << "Pole not initialized!" << std::endl;
       return;
       }
       CopyPole = p->field; //Копия
        
            for (int x = 0; x < (int) p->field[0].size(); ++x) {
                for (int y = 0; y < (int) p->field.size(); ++y) {
                   CopyPole[y][x] = p->field[y][x] > slise ? 255 : 0;
                } 
            }
        
        bmpHandler.bmp_write(CopyPole, filename, logger, b);
        logger.logMessage("Created slised BMP file.", b);
    }
    
    void wave(std::vector<Component>& componenti, std::vector<std::vector<double>> &CopyPole, Pole*& p, int Noise, Logger &logger, bool b) {
    if (p == nullptr) {
       std::cerr << "Pole not initialized!" << std::endl;
       return;
    }
    
            for (unsigned short y = 0; y < (int) p->field.size(); ++y) {
                for (unsigned short x = 0; x < (int) p->field[y].size(); ++x) {
                    if (CopyPole[y][x] <= 255 && CopyPole[y][x] >= 255) {
                        count = 0;
                        Component Componenta(p->field.size(), p->field[y].size());
                        
                        incrementAndCollect(Componenta.componenta, CopyPole, x, y, 1);
                        if (count > Noise) {
                         logger.logMessage("Found component, adding...", b);
                         componenti.emplace_back(Componenta); 
                         logger.logMessage("Component added, amount="  + std::to_string(componenti.size()), b);
                       }
                   }
               }
           }

        logger.logMessage("Wave completed, amount of components = " + std::to_string(componenti.size()), b);
    }
};

class Kmeans {
private:
    double eps = 0.01;
    double h_sum=0;
    double h_j;
    bool boool;
    int indicator; //И ещё один индикатор
    int temp_x, temp_y;
    std::vector<std::vector<float>> center_prev; //для предыдущих центров
    std::vector<std::vector<int>> color; //Карта принадлежности точек кластерам
    int i;
    int i_min;
    double distance;
public:
    
    void k_means(std::vector<std::vector<float>>& center, std::vector<Component>& componenti, Pole*& p, int k, Logger &logger, bool b) 
    {
        
        center.resize(k, std::vector<float>(2, 0));
        center_prev.resize(k, std::vector<float>(2, 0));
        color.resize((int) p->field.size(), std::vector<int>((int) p->field[0].size(), -1));
        
        srand(time(0));  //Шаг 1: выделяем случайные центры
        
        for (i=0;i<k;i++) 
        {                  //Выбираем k случайных центров
            center_prev[i][0] = center[i][0] = rand() % ((int) p->field.size());  
            
            center_prev[i][1] = center[i][1] = rand() % ((int) p->field[0].size()); 
            //std::cout << center[i][0] << ", " << center[i][1] << " from (" << (int) p->field.size() << ", " << (int) p->field[0].size() << ")\n";
        }
                         //Шаг 2: повторяем алгоритм
        
        while (indicator<k)
        {             //подшаг 2.1: определить принадлежность точек кластерам
            
            for (int y = 0; y < (int) p->field.size(); ++y) {
                for (int x = 0; x < (int) p->field[y].size(); ++x) {
                    boool=false; 

                    for (const auto& Componenta : componenti) {if(Componenta.componenta[y][x]>=1) {boool=true; break;}}
                
                    if (boool)
                    {
                        i_min=0;
                        distance=sqrt(pow(y-center[0][0], 2)+pow(x-center[0][1], 2));
                        for (i=0;i<k;i++)
                        {
                            if (sqrt(pow(y-center[i][0], 2)+pow(x-center[i][1], 2))<distance)
                            {
                                distance=sqrt(pow(y-center[i][0], 2)+pow(x-center[i][1], 2));
                                i_min=i;
                            }
                        }
                        color[y][x]=i_min;  
                                                  //std::cout << color[y][x] << " " ;
                        
                    }else { color[y][x]=-1; /*std::cout << "- ";*/ }  
                                                  //if (x==(int) p->field[y].size() - 1) { std::cout<< "" <<std::endl; }
                }
            }         //подшаг 2.2: определить новые центры
            
            for (i=0;i<k;i++)
            {
                h_sum=0;      
                center_prev[i][0]=center[i][0];
                center_prev[i][1]=center[i][1];
                for (int y = 0; y < (int) p->field.size(); ++y) {
                    for (int x = 0; x < (int) p->field[y].size(); ++x) {
                        if (color[y][x]==i)
                        {
                            h_j=0;
                            for (const auto& Componenta : componenti) { h_j+=Componenta.componenta[y][x]; } 
   
                            center[i][0] =(center[i][0]*h_sum+h_j*y)/(h_sum+h_j);
                            center[i][1] =(center[i][1]*h_sum+h_j*x)/(h_sum+h_j); 
                            
                            
                            
                            h_sum+=h_j; 
                        }
                    }
                }
            }
            indicator=0;
            for(i=0;i<k;i++)
            {
                if (sqrt(pow(center_prev[i][0]-center[i][0], 2)+pow(center_prev[i][1]-center[i][1], 2)) < eps ) {indicator++;}
                //std::cout << std::to_string(center[i][0])+" "+std::to_string(center[i][1])+" "+std::to_string(sqrt(pow(center_prev[i][0]-center[i][0], 2)+pow(center_prev[i][1]-center[i][1], 2))) << std::endl;
            }
            
        }
        
        logger.logMessage("Found k centers:", b);
        for (i=0;i<k;i++) {
        logger.logMessage("Center "+std::to_string(i+1)+": ("+std::to_string((int)center[i][1])+", "+std::to_string((int)center[i][0])+")", b); }
    }
    
    
};

class Factors
{
private:
    int len=0;
    int n;
    int k=0;
    double xcp, ycp, xycp, xxcp;
    
    
    
    
public:
    
    void factor_find(std::vector<float>& k1, std::vector<float>& k2, std::vector<float>& b1, std::vector<float>& b2, std::vector<Component>& componenti, Pole*& p, Logger &logger, bool b)
    {
        //for (const auto& Componenta : componenti) {len++;}
        len =(int) componenti.size();
        k1.resize(len, 0);
        k2.resize(len, 0);
        b1.resize(len, 0);
        b2.resize(len, 0);
        
        for (const auto& Componenta : componenti) 
        {
            xcp=ycp=xycp=xxcp=0;
            n=0;
            for (int y = 0; y < (int) p->field.size(); ++y) {
                for (int x = 0; x < (int) p->field[y].size(); ++x) {
                    if (Componenta.componenta[y][x] >=1)
                    {
                        xcp = (xcp*n + x)/(n+1);
                        ycp = (ycp*n + y)/(n+1);
                        xxcp = (xxcp*n + x*x)/(n+1);
                        xycp = (xycp*n + x*y)/(n+1);
                        n++;
                    }
                }
            }
            k1[k]=(xycp-xcp*ycp)/(xxcp-xcp*xcp);
            b1[k]=ycp-k1[k]*xcp;
            k2[k]=-1/k1[k];
            b2[k]=ycp-k2[k]*xcp;
            k++;
        }
        logger.logMessage("All factors found:", b);
        for (int i=0;i<len;i++) {
        logger.logMessage("y="+std::to_string(k1[i])+"*x+"+std::to_string(b1[i])+"; "+"y="+std::to_string(k2[i])+"*x+"+std::to_string(b2[i]), b); }
    }
};

class Triangulator
{
private:
    std::vector<point> points;
    double center_x, center_y;
    double n;
    point p0;
    
    
public:
    void triangulation_build(std::vector<Component>& componenti, std::vector<triangle>& triangles, std::vector<edge>& edges, int width, int height, Logger &logger, bool b)
    {
        points.resize(componenti.size());
        for (int i=0; i<(int)componenti.size(); i++) 
        {
            center_x=0;
            center_y=0;
            n=0;
            for (int y=0; y<(int)componenti[i].componenta.size(); y++)
            {
                for (int x=0; x<(int)componenti[i].componenta[y].size(); x++)
                {
                    if (componenti[i].componenta[y][x]>0)
                    {
                        center_x = (center_x * n + x * componenti[i].componenta[y][x]) / (n+componenti[i].componenta[y][x]);
                        center_y = (center_y * n + y * componenti[i].componenta[y][x]) / (n+componenti[i].componenta[y][x]);
                        n+=componenti[i].componenta[y][x];
                    }
                }
            }
            points[i]=point(center_x, center_y);
            
        }
        
        if (points.size() < 3) {
            logger.logMessage("Недостаточно точек для триангуляции (меньше 3)", b);
            return;
        }

        // Проверка на NaN-точки
        for (const auto& p : points) {
            if (std::isnan(p.x) || std::isnan(p.y)) {
                logger.logMessage("Обнаружена некорректная точка (NaN)", b);
                return;
            }
        }
    
        triangles.emplace_back(point(0, 0), point(0, 2*std::max(width, height)), point(2*std::max(width, height), 0));
        
        
        
        for (const auto& p : points) 
        {
            std::vector<triangle> bad_triangles;
            std::copy_if(triangles.begin(), triangles.end(), std::back_inserter(bad_triangles),
                         [&](const triangle& tri) {                  //точка внутри описанной окружности
                                                      p0 = tri.CircumCenter();
                                                      return ((p.x-p0.x)*(p.x-p0.x)+(p.y-p0.y)*(p.y-p0.y) < (tri.a.x-p0.x)*(tri.a.x-p0.x)+(tri.a.y-p0.y)*(tri.a.y-p0.y));     //расстояние до точки от центра < расстояние от вершины треугольника до центра
                                                });   
            std::vector<edge> polygon;
            
            for (const auto& tri : bad_triangles) {
                for (const auto& edg : { edge{tri.a, tri.b}, edge{tri.b, tri.c}, edge{tri.c, tri.a} }) {
                    bool isShared = std::any_of(bad_triangles.begin(), bad_triangles.end(),
                        [&](const triangle& other) { 
                            return (!(tri == other)) && ((other.a == edg.a && other.b == edg.b) 
                                                    || (other.a == edg.b && other.b == edg.a) 
                                                    || (other.a == edg.a && other.c == edg.b) 
                                                    || (other.a == edg.b && other.c == edg.a) 
                                                    || (other.b == edg.a && other.c == edg.b) 
                                                    || (other.b == edg.b && other.c == edg.a)); });
                    
                    if (!isShared) {
                        polygon.emplace_back(edg);
                        
                    }
                }
            }
            
            triangles.erase(std::remove_if(triangles.begin(), triangles.end(),
                [&](const triangle& t) { 
                    return std::find(bad_triangles.begin(), bad_triangles.end(), t) != bad_triangles.end();
                }), triangles.end());
            
            
            
            for (const auto& edg : polygon) {
                triangles.emplace_back(edg.a, edg.b, p);
            }
            
            
    
        }
        
        triangles.erase(std::remove_if(triangles.begin(), triangles.end(),
            [&](const triangle& t) {
                return t.a == point(0, 0) || t.a == point(0, 2*std::max(width, height)) || t.a == point(2*std::max(width, height), 0) ||
                       t.b == point(0, 0) || t.b == point(0, 2*std::max(width, height)) || t.b == point(2*std::max(width, height), 0) ||
                       t.c == point(0, 0) || t.c == point(0, 2*std::max(width, height)) || t.c == point(2*std::max(width, height), 0); }), triangles.end());
        
        
        
        for (const auto& tri : triangles)
        {
            if (std::find(edges.begin(), edges.end(), edge(tri.a, tri.b)) == edges.end()) { edges.emplace_back(tri.a, tri.b); }
            if (std::find(edges.begin(), edges.end(), edge(tri.b, tri.c)) == edges.end()) { edges.emplace_back(tri.b, tri.c); }
            if (std::find(edges.begin(), edges.end(), edge(tri.a, tri.c)) == edges.end()) { edges.emplace_back(tri.a, tri.c); }
        }
        
        
        std::ofstream triFile("triangulation_output.txt");
        for (const auto& tri : triangles)
        {
            triFile << tri.a.x << " " << tri.a.y << std::endl;
            triFile << tri.b.x << " " << tri.b.y << std::endl;
            triFile << tri.c.x << " " << tri.c.y << "\n" << std::endl;
        }
        triFile.close();
        logger.logMessage("triangulation completed successfully, built " + std::to_string(triangles.size()) + " triangles", b);
        
    }
    
};

class Diagram_builder
{
private:
    std::vector<triangle> tri_stack;
    point cc1, cc2, sr;
    
    
    
public:
    
    void build_diagram(std::vector<Component>& componenti, std::vector<triangle>& triangles, std::vector<edge>& edges, int width, int height, std::vector<edge>&voronoi, Logger &logger, bool b)
    {
        for (const auto& Edge : edges)
        {
            tri_stack.clear();
            for (const auto& tri : triangles)
            {
                if ((tri.a == Edge.a && tri.b == Edge.b) 
                 || (tri.a == Edge.b && tri.b == Edge.a) 
                 || (tri.a == Edge.a && tri.c == Edge.b) 
                 || (tri.a == Edge.b && tri.c == Edge.a) 
                 || (tri.b == Edge.a && tri.c == Edge.b) 
                 || (tri.b == Edge.b && tri.c == Edge.a)) { tri_stack.push_back(tri); }
            }
            cc1=tri_stack[0].CircumCenter();
            if (tri_stack.size()==2) 
            {
                cc2=tri_stack[1].CircumCenter();
                if (cc1.x<0)
                {
                    if (cc1.y-cc1.x*(cc2.y-cc1.y)/(cc2.x-cc1.x)<0) { voronoi.emplace_back(point(cc1.x-cc1.y*(cc2.x-cc1.x)/(cc2.y-cc1.y), 0), cc2); }
                    else if (cc1.y-cc1.x*(cc2.y-cc1.y)/(cc2.x-cc1.x)>height-1) { voronoi.emplace_back(point(cc1.x+(height-1-cc1.y)*(cc2.x-cc1.x)/(cc2.y-cc1.y), height-1), cc2); }
                    else { voronoi.emplace_back(point(0, cc1.y-cc1.x*(cc2.y-cc1.y)/(cc2.x-cc1.x)), cc2); }
                } else if (cc1.y<0)
                {
                    if (cc1.x-cc1.y*(cc2.x-cc1.x)/(cc2.y-cc1.y)>width-1) { voronoi.emplace_back(point(width-1, cc1.y+(width-1-cc1.x)*(cc2.y-cc1.y)/(cc2.x-cc1.x)), cc2); }
                        else { voronoi.emplace_back(point(cc1.x-cc1.y*(cc2.x-cc1.x)/(cc2.y-cc1.y), 0), cc2); }
                } else if (cc1.y>height-1)
                {
                    if (cc1.x+(height-1-cc1.y)*(cc2.x-cc1.x)/(cc2.y-cc1.y)>width-1) { voronoi.emplace_back(point(width-1, cc1.y+(width-1-cc1.x)*(cc2.y-cc1.y)/(cc2.x-cc1.x)), cc2); }
                    else { voronoi.emplace_back(point(cc1.x+(height-1-cc1.y)*(cc2.x-cc1.x)/(cc2.y-cc1.y), height-1), cc2); }
                } else if (cc1.x>width-1)
                {
                    voronoi.emplace_back(point(width-1, cc1.y+(width-1-cc1.x)*(cc2.y-cc1.y)/(cc2.x-cc1.x)), cc2);
                } else if (cc2.x<0)
                {
                    if (cc2.y-cc2.x*(cc1.y-cc2.y)/(cc1.x-cc2.x)<0) { voronoi.emplace_back(point(cc2.x-cc2.y*(cc1.x-cc2.x)/(cc1.y-cc2.y), 0), cc1); }
                    else if (cc2.y-cc2.x*(cc1.y-cc2.y)/(cc1.x-cc2.x)>height-1) { voronoi.emplace_back(point(cc2.x+(height-1-cc2.y)*(cc1.x-cc2.x)/(cc1.y-cc2.y), height-1), cc1); }
                    else { voronoi.emplace_back(point(0, cc2.y-cc2.x*(cc1.y-cc2.y)/(cc1.x-cc2.x)), cc1); }
                } else if (cc2.y<0)
                {
                    if (cc2.x-cc2.y*(cc1.x-cc2.x)/(cc1.y-cc2.y)>width-1) { voronoi.emplace_back(point(width-1, cc2.y+(width-1-cc2.x)*(cc1.y-cc2.y)/(cc1.x-cc2.x)), cc1); }
                    else { voronoi.emplace_back(point(cc2.x-cc2.y*(cc1.x-cc2.x)/(cc1.y-cc2.y), 0), cc1); }
                } else if (cc2.y>height-1)
                {
                    if (cc2.x+(height-1-cc2.y)*(cc1.x-cc2.x)/(cc1.y-cc2.y)>width-1) { voronoi.emplace_back(point(width-1, cc2.y+(width-1-cc2.x)*(cc1.y-cc2.y)/(cc1.x-cc2.x)), cc1); }
                    else { voronoi.emplace_back(point(cc2.x+(height-1-cc2.y)*(cc1.x-cc2.x)/(cc1.y-cc2.y), height-1), cc1); }
                } else if (cc2.x>width-1)
                {
                    voronoi.emplace_back(point(width-1, cc2.y+(width-1-cc2.x)*(cc1.y-cc2.y)/(cc1.x-cc2.x)), cc1);
                } else { voronoi.emplace_back(cc1, cc2); }
            } else
            {
                sr=point((Edge.a.x+Edge.b.x)/2, (Edge.a.y+Edge.b.y)/2);
                if (sr.x-cc1.x>0)
                {
                    if (cc1.y+(width-1-cc1.x)*(sr.y-cc1.y)/(sr.x-cc1.x)<0) { voronoi.emplace_back(point(cc1.x-sr.y*(sr.x-cc1.x)/(sr.y-cc1.y), 0), cc1); }             
                    else if (cc1.y+(width-1-cc1.x)*(sr.y-cc1.y)/(sr.x-cc1.x)>height-1) { voronoi.emplace_back(point(cc1.x+(height-1-sr.y)*(sr.x-cc1.x)/(sr.y-cc1.y), height-1), cc1); }
                    else { voronoi.emplace_back(point(width-1, cc1.y+(width-1-cc1.x)*(sr.y-cc1.y)/(sr.x-cc1.x)), cc1); }
                } else
                {
                    if (cc1.y-cc1.x*(sr.y-cc1.y)/(sr.x-cc1.x)<0) { voronoi.emplace_back(point(cc1.x-sr.y*(sr.x-cc1.x)/(sr.y-cc1.y), 0), cc1); }             
                    else if (cc1.y-cc1.x*(sr.y-cc1.y)/(sr.x-cc1.x)>height-1) { voronoi.emplace_back(point(cc1.x+(height-1-sr.y)*(sr.x-cc1.x)/(sr.y-cc1.y), height-1), cc1); }
                    else { voronoi.emplace_back(point(0, cc1.y-cc1.x*(sr.y-cc1.y)/(sr.x-cc1.x)), cc1); }
                }
            }
            
        }
        
        std::ofstream voronoiFile("voronoi_output.txt");
        for (const auto& edg : voronoi)
        {
            voronoiFile << edg.a.x << " " << edg.a.y << std::endl;
            voronoiFile << edg.b.x << " " << edg.b.y << "\n" << std::endl;
        }
        voronoiFile.close();
        
        logger.logMessage("diagram built successfully", b);
        
    }
    
};

class node
{
public:
    point p;
    double g;
    double h;
    bool visited;
    node* prev;
    
    node(point p_ = point(), double g_ = 100000000, double h_ = 100000000, bool visited_ = false): p(p_), g(g_), h(h_), visited(visited_) {prev = nullptr;}
    
};

class Pathfinder
{
private:
    
    double k, b, f;
    double x1, x2, y1, y2;
    bool radius_failure;
    double max_roll, max_pitch;
    point start_p, finish_p;
    double min_dist_start, min_dist_finish;
    std::vector<node> graph;
    
    point temp_point;
    
    node* curr;
    double new_g;
    
public:
    
    std::vector<node>::iterator find_node(std::vector<node>::iterator first, std::vector<node>::iterator last, point pnt)
    {
        for (; first != last; ++first)
        {
            if ((*first).p == pnt) { return first; }
        }
        return last;
    }
    
    
    void pathfind(std::vector<Component>& componenti, Pole*& p, int slise, std::vector<edge>&voronoi, std::vector<edge>&Path, double x0, double y0, double xk, double yk, double R, double pitch, double roll, Logger &logger, bool b)
    {
        std::vector<edge> edges;
        std::copy(voronoi.begin(), voronoi.end(), std::back_inserter(edges));
        if (x0<0 || x0>p->field[0].size()-1 || y0<0 || y0>-1+p->field.size())
        {
            logger.logMessage("Incorrect starting point!", b);
            return;
        }
        if (xk<0 || xk>p->field[0].size()-1 || yk<0 || yk>-1+p->field.size())
        {
            logger.logMessage("Incorrect finish point!", b);
            return;
        }
        for (int y = 0; y < (int) p->field.size(); ++y) 
        {
            for (int x = 0; x < (int) p->field[y].size(); ++x) 
            {
                if ((p->field[y][x]>slise) && ((x-x0)*(x-x0)+(y-y0)*(y-y0)<R*R))
                {
                    logger.logMessage("Incorrect starting point!", b);
                    return;
                }
                if ((p->field[y][x]>slise) && ((x-xk)*(x-xk)+(y-yk)*(y-yk)<R*R))
                {
                    logger.logMessage("Incorrect finish point!", b);
                    return;
                }
            }
        }
        for(auto edgi{edges.end()-1};;)
        {
            radius_failure=false;
            max_roll=0;
            max_pitch=0;
            k=((*edgi).b.y-(*edgi).a.y)/((*edgi).b.x-(*edgi).a.x);
            b=(*edgi).b.y-k*((*edgi).b.x);
            for (int y = 0; y < (int) p->field.size(); ++y) 
            {
                for (int x = 0; x < (int) p->field[y].size(); ++x) 
                {
                    f=k*x-y+b;
                    if ( (p->field[y][x]>slise) && (f<R*sqrt(k*k+1) || f>-R*sqrt(k*k+1)) && ((x-f/(2*k)>(*edgi).a.x && x-f/(2*k)<(*edgi).b.x) || (x-f/(2*k)<(*edgi).a.x && x-f/(2*k)>(*edgi).b.x)) && ((y+f/2>(*edgi).a.y && y+f/2<(*edgi).b.y) || (y+f/2<(*edgi).a.y && y+f/2>(*edgi).b.y)) )   //проверка на застревание при прохождении по ребру
                    {
                        
                        radius_failure=true;
                        break;
                    }
                    if (int(y+f/2)>0 && int(y+f/2)<(int) p->field.size()-1 && int(x-f/(2*k))>0 && int(x-f/(2*k))< (int) p->field[y].size()-1)
                    {
                        if ((max_roll<std::abs(atan((p->field[y][x]-p->field[int(y+f/2)][int(x-f/(2*k))])/((f/(2*k))*(f/(2*k))+(f/2)*(f/2))))) && (f<R*sqrt(k*k+1) || f>-R*sqrt(k*k+1)) && ((x-f/(2*k)>(*edgi).a.x && x-f/(2*k)<(*edgi).b.x) || (x-f/(2*k)<(*edgi).a.x && x-f/(2*k)>(*edgi).b.x)) && ((y+f/2>(*edgi).a.y && y+f/2<(*edgi).b.y) || (y+f/2<(*edgi).a.y && y+f/2>(*edgi).b.y)))   //проверка на боковой угол
                        {
                            max_roll=std::abs(atan( (p->field[y][x] - p->field[int(y+f/2)][int(x-f/(2*k))]) / ((f/(2*k))*(f/(2*k))+(f/2)*(f/2)) ));
                        }
                    }
                }
                if (radius_failure) {break;}
                if (y>R && y<p->field.size()-R-1 && (y>(*edgi).a.y && y<(*edgi).b.y) || (y<(*edgi).a.y && y>(*edgi).b.y)) 
                {
                    y1=y-R;
                    y2=y+R;
                    x1=(y1-b)/k;
                    x2=(y2-b)/k;
                    if (int(x1)>0 && int(x1)<p->field[y].size()-1 && int(x2)>0 && int(x2)<p->field[y].size()-1)
                    {
                        if (max_pitch < std::abs(atan((p->field[int(y1)][int(x1)]-p->field[int(y2)][int(x2)])/((y1-y2)*(y1-y2)+(x1-x2)*(x1-x2)))) )
                        {
                            max_pitch = std::abs(atan((p->field[int(y1)][int(x1)]-p->field[int(y2)][int(x2)])/((y1-y2)*(y1-y2)+(x1-x2)*(x1-x2))));
                        }
                    }
                }
            }
            if (radius_failure) 
            {
                if (edgi == edges.begin()) {edges.erase(edgi); logger.logMessage("Edge (("+std::to_string((*edgi).a.x)+", "+std::to_string((*edgi).a.y)+"), ("+std::to_string((*edgi).b.x)+"), ("+std::to_string((*edgi).b.y)+")) doesnt fit radius requirments", b); break;}
                else
                {
                --edgi;
                edges.erase(edgi+1);
                }
                logger.logMessage("Edge (("+std::to_string((*edgi).a.x)+", "+std::to_string((*edgi).a.y)+"), ("+std::to_string((*edgi).b.x)+"), ("+std::to_string((*edgi).b.y)+")) doesnt fit radius requirments", b);
                continue;
            } else if (max_roll*180/3.141592653>roll) 
            {
                if (edgi == edges.begin()) {edges.erase(edgi); logger.logMessage("Edge (("+std::to_string((*edgi).a.x)+", "+std::to_string((*edgi).a.y)+"), ("+std::to_string((*edgi).b.x)+"), ("+std::to_string((*edgi).b.y)+")) doesnt fit roll requirments", b); break;}
                else
                {
                --edgi;
                edges.erase(edgi+1);
                }
                logger.logMessage("Edge (("+std::to_string((*edgi).a.x)+", "+std::to_string((*edgi).a.y)+"), ("+std::to_string((*edgi).b.x)+"), ("+std::to_string((*edgi).b.y)+")) doesnt fit roll requirments", b);
                continue;
            } else if (max_pitch*180/3.141592653>pitch) 
            {
                if (edgi == edges.begin()) {edges.erase(edgi); logger.logMessage("Edge (("+std::to_string((*edgi).a.x)+", "+std::to_string((*edgi).a.y)+"), ("+std::to_string((*edgi).b.x)+"), ("+std::to_string((*edgi).b.y)+")) doesnt fit pitch requirments", b); break;}
                else
                {
                --edgi;
                edges.erase(edgi+1);
                }
                logger.logMessage("Edge (("+std::to_string((*edgi).a.x)+", "+std::to_string((*edgi).a.y)+"), ("+std::to_string((*edgi).b.x)+"), ("+std::to_string((*edgi).b.y)+")) doesnt fit pitch requirments", b);
                continue;
            }
            if (edgi == edges.begin()) {break;}
            else { --edgi; }
        }
        min_dist_start = p->field.size() * p->field[0].size();
        min_dist_finish = p->field.size() * p->field[0].size();
        start_p = edges[0].a;
        finish_p = edges[0].a;
        for(auto edgi{edges.begin()}; edgi !=edges.end(); ++edgi )
        {
            if (min_dist_start > sqrt(((*edgi).a.x-x0)*((*edgi).a.x-x0)+((*edgi).a.y-y0)*((*edgi).a.y-y0)))
            {
                min_dist_start=sqrt(((*edgi).a.x-x0)*((*edgi).a.x-x0)+((*edgi).a.y-y0)*((*edgi).a.y-y0));
                start_p = (*edgi).a;
            }
            if (min_dist_start > sqrt(((*edgi).b.x-x0)*((*edgi).b.x-x0)+((*edgi).b.y-y0)*((*edgi).b.y-y0)))
            {
                min_dist_start=sqrt(((*edgi).b.x-x0)*((*edgi).b.x-x0)+((*edgi).b.y-y0)*((*edgi).b.y-y0));
                start_p = (*edgi).b;
            }
            if (min_dist_finish > sqrt(((*edgi).a.x-xk)*((*edgi).a.x-xk)+((*edgi).a.y-yk)*((*edgi).a.y-yk)))
            {
                min_dist_finish=sqrt(((*edgi).a.x-xk)*((*edgi).a.x-xk)+((*edgi).a.y-yk)*((*edgi).a.y-yk));
                finish_p = (*edgi).a;
            } 
            if (min_dist_finish > sqrt(((*edgi).b.x-xk)*((*edgi).b.x-xk)+((*edgi).b.y-yk)*((*edgi).b.y-yk)))
            {
                min_dist_finish=sqrt(((*edgi).b.x-xk)*((*edgi).b.x-xk)+((*edgi).b.y-yk)*((*edgi).b.y-yk));
                finish_p = (*edgi).b;
            }
            
        }
        for(const auto& edg : {edge{start_p, point{x0, y0}}, edge{finish_p, point{xk, yk}}})      //проверка новых рёбер
        {
            radius_failure=false;
            max_roll=0;
            max_pitch=0;
            k=(edg.b.y-edg.a.y)/(edg.b.x-edg.a.x);
            b=edg.b.y-k*(edg.b.x);
            for (int y = 0; y < (int) p->field.size(); ++y) 
            {
                for (int x = 0; x < (int) p->field[y].size(); ++x) 
                {
                    f=k*x-y+b;
                    if ( (p->field[y][x]>slise) && (f<R*sqrt(k*k+1) || f>-R*sqrt(k*k+1)) && ((x-f/(2*k)>edg.a.x && x-f/(2*k)<edg.b.x) || (x-f/(2*k)<edg.a.x && x-f/(2*k)>edg.b.x)) && ((y+f/2>edg.a.y && y+f/2<edg.b.y) || (y+f/2<edg.a.y && y+f/2>edg.b.y)) )   //проверка на застревание при прохождении по ребру
                    {
                        
                        radius_failure=true;
                        break;
                    }
                    if (int(y+f/2)>0 && int(y+f/2)<(int) p->field.size()-1 && int(x-f/(2*k))>0 && int(x-f/(2*k))< (int) p->field[y].size()-1)
                    {
                        if ((max_roll<std::abs(atan((p->field[y][x]-p->field[int(y+f/2)][int(x-f/(2*k))])/((f/(2*k))*(f/(2*k))+(f/2)*(f/2))))) && (f<R*sqrt(k*k+1) || f>-R*sqrt(k*k+1)) && ((x-f/(2*k)>edg.a.x && x-f/(2*k)<edg.b.x) || (x-f/(2*k)<edg.a.x && x-f/(2*k)>edg.b.x)) && ((y+f/2>edg.a.y && y+f/2<edg.b.y) || (y+f/2<edg.a.y && y+f/2>edg.b.y)))   //проверка на боковой угол
                        {
                            max_roll=std::abs(atan( (p->field[y][x] - p->field[int(y+f/2)][int(x-f/(2*k))]) / ((f/(2*k))*(f/(2*k))+(f/2)*(f/2)) ));
                        }
                    }
                }
                if (radius_failure) {break;}
                if (y>R && y<p->field.size()-R-1 && ((y>edg.a.y && y<edg.b.y) || (y<edg.a.y && y>edg.b.y))) 
                {
                    y1=y-R;
                    y2=y+R;
                    x1=(y1-b)/k;
                    x2=(y2-b)/k;
                    if (int(x1)>0 && int(x1)<p->field[y].size()-1 && int(x2)>0 && int(x2)<p->field[y].size()-1)
                    {
                        if ((max_pitch < std::abs(atan((p->field[int(y1)][int(x1)]-p->field[int(y2)][int(x2)])/((y1-y2)*(y1-y2)+(x1-x2)*(x1-x2))))) )
                        {
                            max_pitch = std::abs(atan((p->field[int(y1)][int(x1)]-p->field[int(y2)][int(x2)])/((y1-y2)*(y1-y2)+(x1-x2)*(x1-x2))));
                        }
                    }
                }
            }
            if (radius_failure) 
            {
                logger.logMessage("One of added edges doesnt fit radius requirments", b);
                return;
            } else if (max_roll*180/3.141592653>roll) 
            {
                logger.logMessage("One of added edges doesnt fit roll requirments", b);
                return;
            } else if (max_pitch*180/3.141592653>pitch) 
            {
                logger.logMessage("One of added edges doesnt fit pitch requirments", b);
                return;
            }
        }
        Path.emplace_back(point{x0, y0}, start_p);                 //добавляем стартовое ребро
        for (const auto& edg : edges)                                               //строим узлы по рёбрам
        {
            if (find_node(graph.begin(), graph.end(), edg.a) == graph.end()) {graph.emplace_back(edg.a);}
            if (find_node(graph.begin(), graph.end(), edg.b) == graph.end()) {graph.emplace_back(edg.b);}
        }                                                                       //устанавливаем значения g и h для начальной точки
        (*find_node(graph.begin(), graph.end(), start_p)).g=0;
        (*find_node(graph.begin(), graph.end(), start_p)).h=0;
        
        auto current_node { graph.begin() };                                   // тип - итератор!
        while (true)
        {
            
            current_node=graph.begin();
            for (auto nodi{graph.begin()}; nodi != graph.end(); nodi++ )                         //ищем узел с наименьшей f
            {
                if ((*nodi).visited==false && (*nodi).g+(*nodi).h<(*current_node).g+(*current_node).h)
                {
                    current_node=nodi;
                }
            }
            
            if ((*current_node).p==finish_p)                                                 //если конечная точка - завершаем программу
            {
                curr = &*current_node;
                while (curr->prev!=nullptr)
                {
                    Path.emplace_back(curr->prev->p, curr->p);
                    curr = curr->prev;
                }
                break;
            } 
            (*current_node).visited = true;
            
            for (const auto& edg : edges)                     // обновляем g и h для соседей
            {
                if (edg.a == (*current_node).p) 
                {
                    auto next_node{find_node(graph.begin(), graph.end(), edg.b)};
                    new_g=(*current_node).g + sqrt(  ((*current_node).p.x-(*next_node).p.x)*((*current_node).p.x-(*next_node).p.x)+((*current_node).p.y-(*next_node).p.y)*((*current_node).p.y-(*next_node).p.y)  );
                    if (new_g < (*next_node).g)
                    {
                        (*next_node).g=new_g;
                        (*next_node).h=sqrt(  (finish_p.x-(*next_node).p.x)*(finish_p.x-(*next_node).p.x)+(finish_p.y-(*next_node).p.y)*(finish_p.y-(*next_node).p.y)  );
                        (*next_node).prev = &*current_node;
                    }
                    
                    
                }
                if (edg.b == (*current_node).p) 
                {
                    auto next_node{find_node(graph.begin(), graph.end(), edg.a)};
                    new_g=(*current_node).g + sqrt(  ((*current_node).p.x-(*next_node).p.x)*((*current_node).p.x-(*next_node).p.x)+((*current_node).p.y-(*next_node).p.y)*((*current_node).p.y-(*next_node).p.y)  );
                    if (new_g < (*next_node).g)
                    {
                        (*next_node).g=new_g;
                        (*next_node).h=sqrt(  (finish_p.x-(*next_node).p.x)*(finish_p.x-(*next_node).p.x)+(finish_p.y-(*next_node).p.y)*(finish_p.y-(*next_node).p.y)  );
                        (*next_node).prev = &*current_node;
                    }
                    
                    
                }
            }
            
            current_node=graph.begin();
            for (auto nodi{graph.begin()}; nodi != graph.end(); nodi++ )
            {
                if ((*nodi).visited==false && (*nodi).g+(*nodi).h<(*current_node).g+(*current_node).h)
                {
                    current_node=nodi;
                    
                }
            }
            if ((*current_node).g > 99999999) 
            {
                logger.logMessage("failed to find the path", b);
                return;
            }
        } 
        
        Path.emplace_back(finish_p, point{xk, yk});
        
        std::ofstream PathFile("path_output.txt");
        for (const auto& edg : Path)
        {
            PathFile << edg.a.x << " " << edg.a.y << std::endl;
            PathFile << edg.b.x << " " << edg.b.y << "\n" << std::endl;
        }
        PathFile.close();
        
        logger.logMessage("Path built successfully", b);
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
    Kmeans kmeans;
    Factors factors;
    Triangulator triangulator;
    Diagram_builder diagram_builder;
    Pathfinder pathfinder;
    
    Pole* p = nullptr; // Pointer to Pole
    std::vector<std::vector<float>> centers;  //массив центров c 2 координатами, y и x
    std::vector<float> k1;
    std::vector<float> k2;
    std::vector<float> b1;
    std::vector<float> b2;
    std::vector<triangle> triangles;
    std::vector<edge> edges;
    std::vector<edge> voronoi;
    std::vector<edge> Path;
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
    
void Dispetcher(const std::string& command, int A, int B, double h, double x, double y, double sigma_x, double sigma_y, 
std::vector<std::vector<double>>& pixelMatrix, const std::string& filename, int slise, int k, int Noise, double R, double pitch, double roll, double x0, double y0, double xk, double yk) {
   
    if (command == "init") {
    gaussBuilder.init(A, B, p, loggercontrol, b);
    }

    if (command == "g") {
    gaussBuilder.addgauss(h, x, y, sigma_x, sigma_y, gaussi, loggercontrol, b);
    }
    
    if (command == "generate") {
    gaussBuilder.generate(p, gaussi, loggercontrol, b);
    }
    
    if (command == "gnuplot") {
    gnuplotInterface.gnuplot(p, loggercontrol, b);
    }

    if (command == "bmp_write") {
    bmpHandler.bmp_write(pixelMatrix, filename, loggercontrol, b);
    }

    if (command == "bmp_read") {
    bmpHandler.bmp_read(gaussBuilder, filename, pixelMatrix, p, loggercontrol, b);
    }

    if (command == "bin") {
    componentCalculator.bin(bmpHandler, CopyPole, slise, filename, p, loggercontrol, b);
    }
    
    if (command == "wave") {
    componentCalculator.wave(componenti, CopyPole, p, Noise, loggercontrol, b);
    }
    
    if (command == "k_means") {
    kmeans.k_means(centers, componenti, p, k, loggercontrol, b);
    }
    
    if (command == "find_factors") {
    factors.factor_find(k1, k2, b1, b2, componenti, p, loggercontrol, b);
    }
    
    if (command == "triangulate") {
    triangulator.triangulation_build(componenti, triangles, edges, A, B, loggercontrol, b);
    }
    
    if (command == "voronoi") {
    diagram_builder.build_diagram(componenti, triangles, edges, A, B, voronoi, loggercontrol, b);
    }
    if (command == "pathfind") {
    pathfinder.pathfind(componenti, p, slise, voronoi, Path, x0, y0, xk, yk, R, pitch, roll, loggercontrol, b);
    }
}
};
  


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
        double x = 0, y = 0, sx = 0, sy = 0, h = 0, x0 = 0, y0 = 0, xk = 0, yk = 0;
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
                loggerinterface.logMessage("Error: File not found.", config.loggingInterfaceEnabled);
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
                helpFile << "g [x] [y] [sigma_x] [sigma_y] [h] - создает гаусс с 5 параметрами (int, int, double, double, double)\n";
                helpFile << "generate - складывает гауссы\n";
                helpFile << "gnuplot - рисует картинку в gnuplot\n";
                helpFile << "bmp_write [outputfilename.bmp]- создает черно-белую серую картинку bmp\n";
                helpFile << "bmp_read [inputfilename.bmp]- чтение bmp файла и инициализация поля новыми размерами\n";
                helpFile << "bin [slise] [outputfilename.bmp]- срез: все что выше черное, все что ниже белое\n";
                helpFile << "wave - считается число компонент, каждая компонента запомнена\n";
                helpFile << "k_means [k] - выделение k кластеров\n";
                helpFile << "find_factors - считаются главные и побочные факторы всех компонент\n";
                helpFile << "triangulate - триангуляция на центрах компонент\n";
                helpFile << "voronoi - построение диаграмы Вороного (необходимо перед этим написать triangulate!)\n";
                helpFile << "pathfind [x0] [y0] [xk] [yk] - прокладывает путь из точки (x0, y0) в точку (xk, yk). Для работы необходимо перед этим задать срез через bin !!!\n";
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
    c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, "output.bmp", slice, k, config.Noise, config.R, config.pitch, config.roll, x0, y0, xk, yk);
    loggerinterface.logMessage("Field initialized.", b);
} else if (n != 1) {
                     std::cout << "The init command was not used.\nError\n";
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
                c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, "output.bmp", slice, k, config.Noise, config.R, config.pitch, config.roll, x0, y0, xk, yk);
            } else if (s == "generate") {
                    c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, "output.bmp", slice, k, config.Noise, config.R, config.pitch, config.roll, x0, y0, xk, yk);
                    loggerinterface.logMessage("Generated values in the field.", b);
                } else if (s == "gnuplot") {
                    c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, "output.bmp", slice, k, config.Noise, config.R, config.pitch, config.roll, x0, y0, xk, yk);
                    loggerinterface.logMessage("Called gnuplot.", b);
                } else if (s == "bmp_write") {
                    file >> filename; // Чтение имени файла для bmp_write
                    c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, filename, slice, k, config.Noise, config.R, config.pitch, config.roll, x0, y0, xk, yk);
                    loggerinterface.logMessage("Created BMP file: " + filename, b);
                } else if (s == "bmp_read") {
                    file >> filename; // Чтение имени файла для bmp_read
                    c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, filename, slice, k, config.Noise, config.R, config.pitch, config.roll, x0, y0, xk, yk);
                    loggerinterface.logMessage("Read BMP file: " + filename, b);
                } else if (s == "bin") {
                    file >> slice;
                    file >> filename; // Чтение имени файла для slice
                    c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, filename, slice, k, config.Noise, config.R, config.pitch, config.roll, x0, y0, xk, yk);
                    loggerinterface.logMessage("Slice applied: slice=" + std::to_string(slice), b);
                } else if (s == "wave") {
                    loggerinterface.logMessage("Wave will be used", b);
                    c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, filename, slice, k, config.Noise, config.R, config.pitch, config.roll, x0, y0, xk, yk);
                } else if (s == "k_means") {
                    file >> k; 
                    c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, filename, slice, k, config.Noise, config.R, config.pitch, config.roll, x0, y0, xk, yk);
                    loggerinterface.logMessage("Successfully found " + std::to_string(k) + " means", b);
                }  else if (s == "find_factors") {
                    loggerinterface.logMessage("Finding factors...", b);
                    c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, filename, slice, k, config.Noise, config.R, config.pitch, config.roll, x0, y0, xk, yk);
                }  else if (s == "triangulate") {
                    loggerinterface.logMessage("Triangulating...", b);
                    c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, filename, slice, k, config.Noise, config.R, config.pitch, config.roll, x0, y0, xk, yk);
                }  else if (s == "voronoi") {
                    loggerinterface.logMessage("Started building voronoi diagram", b);
                    c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, filename, slice, k, config.Noise, config.R, config.pitch, config.roll, x0, y0, xk, yk);
                }  else if (s == "pathfind") {
                    file >> x0;
                    file >> y0;
                    file >> xk;
                    file >> yk;
                    loggerinterface.logMessage("Pathfinding started", b);
                    c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, filename, slice, k, config.Noise, config.R, config.pitch, config.roll, x0, y0, xk, yk);
                }
            }
        } else {
            int n = 0;
            while (true) {
                std::cout << "Enter command (help, init, g, generate, gnuplot, bmp_write, bmp_read, bin, wave, k_means, find_factors, triangulate, voronoi, pathfind, end)->";
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
                helpFile << "g [x] [y] [sigma_x] [sigma_y] [h] - создает гаусс с 5 параметрами (int, int, double, double, double)\n";
                helpFile << "generate - складывает гауссы\n";
                helpFile << "gnuplot - рисует картинку в gnuplot\n";
                helpFile << "bmp_write [outputfilename.bmp]- создает черно-белую серую картинку bmp\n";
                helpFile << "bmp_read [inputfilename.bmp]- чтение bmp файла и инициализация поля новыми размерами\n";
                helpFile << "bin [slise] [outputfilename.bmp]- срез: все что выше черное, все что ниже белое\n";
                helpFile << "wave - считается число компонент, каждая компонента запомнена\n";
                helpFile << "k_means [k] - выделение k кластеров\n";
                helpFile << "find_factors - считаются главные и побочные факторы всех компонент\n";
                helpFile << "triangulate - триангуляция на центрах компонент\n";
                helpFile << "voronoi - построение диаграмы Вороного (необходимо перед этим написать triangulate!)\n";
                helpFile << "pathfind [x0] [y0] [xk] [yk] - прокладывает путь из точки (x0, y0) в точку (xk, yk). Для работы необходимо перед этим задать срез через bin !!!\n";
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
    c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, "output.bmp", slice, k, config.Noise, config.R, config.pitch, config.roll, x0, y0, xk, yk);
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
               c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, "output.bmp", slice, k, config.Noise, config.R, config.pitch, config.roll, x0, y0, xk, yk);
            }

                  if (s == "generate") {
                    c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, "output.bmp", slice, k, config.Noise, config.R, config.pitch, config.roll, x0, y0, xk, yk);
                    std::cout << "Generated values in the field." << std::endl;
                    loggerinterface.logMessage("Generated values in the field.", b);
                } if (s == "gnuplot") {
                    c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, "output.bmp", slice, k, config.Noise, config.R, config.pitch, config.roll, x0, y0, xk, yk);
                    std::cout << "Called gnuplot." << std::endl;
                    loggerinterface.logMessage("Called gnuplot.", b);
                } if (s == "bmp_write") {
                    std::cout << "Enter the output filename:" << std::endl;
                    std::cin >> filename;
                    c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, filename, slice, k, config.Noise, config.R, config.pitch, config.roll, x0, y0, xk, yk);
                    std::cout << "Created BMP file." << std::endl;
                    loggerinterface.logMessage("Created BMP file.", b);
                } if (s == "bmp_read") {
                    std::cout << "Enter the filename to read:" << std::endl;
                    std::cin >> filename;
                    c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, filename, slice, k, config.Noise, config.R, config.pitch, config.roll, x0, y0, xk, yk);
                    std::cout << "Read BMP file: " + filename << std::endl;
                    loggerinterface.logMessage("Read BMP file: " + filename, b);
                }  if (s == "end") {
                    std::cout << "Ending the program" << std::endl;
                    loggerinterface.logMessage("Ending the program.", b);
                    break;
                }
                   if (s == "bin") {
                     std::cout << "Enter slice level->" << std::endl;
                     std::cin >> slice;
                     std::cout << "Enter the output filename->" << std::endl;
                     std::cin >> filename;
                     c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, filename, slice, k, config.Noise, config.R, config.pitch, config.roll, x0, y0, xk, yk);
                     loggerinterface.logMessage("Slice applied: slice=" + std::to_string(slice), b);
                     std::cout << "Slice applied: slice=" << slice << std::endl;
                }
                   if (s == "wave") {
                     loggerinterface.logMessage("Wave will be used", b);
                     c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, filename, slice, k, config.Noise, config.R, config.pitch, config.roll, x0, y0, xk, yk);
                }
                   if (s == "k_means") {
                     std::cout << "Enter amount of clusters for k_means->" << std::endl;
                     std::cin >> k;
                     c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, filename, slice, k, config.Noise, config.R, config.pitch, config.roll, x0, y0, xk, yk);
                     loggerinterface.logMessage("Successfully found " + std::to_string(k) + " means", b); 
                }  
                   if (s == "find_factors") {
                     loggerinterface.logMessage("Finding factors...", b);
                     c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, filename, slice, k, config.Noise, config.R, config.pitch, config.roll, x0, y0, xk, yk);
                }  if (s == "triangulate") {
                     loggerinterface.logMessage("Triangulating...", b);
                     c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, filename, slice, k, config.Noise, config.R, config.pitch, config.roll, x0, y0, xk, yk);
                }  if (s == "voronoi") {
                     loggerinterface.logMessage("Started building voronoi diagram", b);
                     c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, filename, slice, k, config.Noise, config.R, config.pitch, config.roll, x0, y0, xk, yk);
                }  if (s == "pathfind") {
                     std::cout<< "Enter starting x and y and final x and y divided by spaces ->" << std::endl;
                     std::cin >> x0 >> y0 >> xk >> yk;
                     loggerinterface.logMessage("Pathfinding started", b);
                     c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, filename, slice, k, config.Noise, config.R, config.pitch, config.roll, x0, y0, xk, yk);
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
