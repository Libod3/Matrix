#include <iostream>
#include <vector>
#include <string>
#include <cmath>
//#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
//#include "doctest.h"


class Matrix { // - Тело класса
private:
    //Вид нашей матрицы
    std::vector<std::vector<double>> data;
public:
//_____________________________________________________________________________________________________________________
    //Просто так создаём пустую матрицу
    Matrix() {}

    // Конструктор матриц, где: row - длина, coloms - ширина, value_of_element - значения элементов
    Matrix(int rows, int cols, double value_of_element) {
        data.resize(rows);
        for (auto &row: data) {
            row.resize(cols, value_of_element);
        }
    }

    // Конструктор матриц, которые на вход получают вектор векторов double
    Matrix(const std::vector<std::vector<double>> &matrix_data) : data(matrix_data) {}

    //ПО ТЗ функция копирования матрицы.
    Matrix(const Matrix &other) : data(other.data) {}
    // Сигнатура конструктора:
    //Matrix(const Matrix& other): Конструктор принимает константную ссылку на другой объект Matrix. Использование
    // константной ссылки важно, так как это предотвращает модификацию объекта-аргумента и обеспечивает
    // эффективность, поскольку объект не копируется полностью в память при передаче в конструктор.
    //Инициализация членов:
    //: data(other.data): Это список инициализации членов класса. Здесь член data нового объекта инициализируется
    // копией данных data из объекта other. Так как data является std::vector<std::vector<double>>, эта операция
    // копирует все векторы и их содержимое. Векторы в C++ управляют своей памятью, поэтому когда вы копируете
    // вектор, создаются новые копии всех содержащихся в нём элементов, а не просто копия ссылок или указателей.
//_____________________________________________________________________________________________________________________

    //size_t - специальный тип данных, который в C++ предназначен для вывода столбцов и индексов
    // для начала буду использовать int
    int get_row() { // функция для получения строк, отлично подойдет .size()
        return data.size();
    }

    int get_coloms() { //Функция вывода столбцов
        //Важно помнить, что если матрица пустая, то вызов нулевого элемента закончится плохо.
        //Делаем проверку на непустоту матрицы. Используем empty()
        if (data.empty()) {
            return 0;
        }
        return data[0].size();
    }

//_____________________________________________________________________________________________________________________
    //Теперь необходимо сделать функцию вывода матрицы.
    void print() {
        // Чисто косметика чтобы не сливались матрицы
        for (int i = 0; i < (2 * data[0].size() - 1); ++i) {
            std::cout << "-";
        }
        std::cout << std::endl;

        for (const auto &row: data) {
            for (double elem: row) {
                std::cout << elem << " ";
            }
            std::cout << std::endl;
        }
        // Тоже чисто для красоты
        for (int i = 0; i < (2 * data[0].size() - 1); ++i) {
            std::cout << "-";
        }
        std::cout << std::endl;
    }

//_____________________________________________________________________________________________________________________
//  Чисто чтоб рисовать матрицами
    void print_cat() const {
        for (const auto& row : data) {
            for (int elem : row) {
                std::cout << (elem ? '*' : ' ') << " ";
            }
            std::cout << std::endl;
        }
    }

//_____________________________________________________________________________________________________________________
    // Функции для нахождения определителя матрицы. Я выбрал именно это способ только из-за того, что использовал
    // именно его про подсчёте определителя "НА БУМАГЕ". Отдельное спасибо Ворощук Анне за знакомство с ним.
    //  В интернете зачастую реализуют метод Гаусса.

    // функция для подсчёта определителя
    double determinant() const {
        if (data.size() != data[0].size()){
            throw std::invalid_argument("Матрицы не квадратная");
            } else {
            int n = data.size(); // Тк матрица квадратная, нам достаточно знать только кол-во строк
            if (n == 1) {
                return data[0][0];
                }
            if (n == 2) {
                return data[0][0] * data[1][1] - data[0][1] * data[1][0];
                }

            double det = 0;
            for (int j = 0; j < n; j++) {
                det += data[0][j] * cofactor(0, j);
                }
            return det;
            }

        }

    // Функция для вычисления алгебраического дополнения
    double cofactor(int row, int coloms) const {
        return std::pow(-1, row + coloms) * minor(row, coloms).determinant();
    }
    //pow нужна для возведения числа в степень (для чередования знака при подсчёте определителя матрицы)
    //#include <cmath> именно для std::pow

    //Функция для вычисления миноров
    Matrix minor(int row, int coloms) const {
        std::vector<std::vector<double>> minorData;
        int n = data.size();
        minorData.reserve(n - 1); // Сделано для оптимизации (установка минимального размера вектора (строк в матрице))

        for (int i = 0; i < n; i++) {
            if (i == row) continue; // continue тут используется для того, чтобы пропустить строку и создать минор
            std::vector<double> row_data;
            row_data.reserve(n - 1);
            for (int j = 0; j < n; j++) {
                if (j == coloms) continue; // continue тут используется для того, чтобы пропустить столбец и создать
                                           // минор (кек копипаст коментария ;З)
                row_data.push_back(data[i][j]);
            }
            minorData.push_back(row_data);
        }
        return Matrix(minorData);
    }
//_____________________________________________________________________________________________________________________
    // Тут я реализую функцию поиска элемента по индексам
    double find_elem(int i, int j){
        if (i >= data.size()){
            std::cout << "Вы ввели слишком большой индекс строки" << " ";
            return 0;
        }
        else if (j >= data[0].size()){
            std::cout << "Вы ввели слишком большой индекс столбца" << " ";
            return 0;
        }
        else{
            return data[i][j];
        }
    }
//_____________________________________________________________________________________________________________________
    // Функция замены элемента в матрице
    void set_elem(int i, int j, double num){
        if (i < data.size() && i >= 0){
            if (j < data[0].size() && j >= 0){
                data[i][j] = num;
            } else {
                throw std::invalid_argument("Такого столбца нет в матрице");
            }
        } else {
            throw std::invalid_argument("Такой строки нет в матрице");
        }
    }
//_____________________________________________________________________________________________________________________
    //Эта функция будет отвечать за сложение двух матриц.
    Matrix sum_mtx(const Matrix& other) const {
        // собираем данные для того чтобы не вызывать функцию повторно
        int row1 = data.size();
        int row2 = other.data.size();
        int coloms1 = data[0].size();
        int coloms2 = other.data[0].size();
        // Также создаём два вектора, один из которых мы будем использовать как строки другого,
        // а второй будет результатом выполнения функции
        std::vector<std::vector<double>> result;
        std::vector<double> result_row;


        // Переходим к проверке условий и реализации
        if (row1 == row2){
            if (coloms1 == coloms2){
                for (int i = 0; i < row1; i++){
                    for (int j = 0; j < coloms1; j++){
                        double elem = data[i][j] + other.data[i][j];
                        result_row.push_back(elem);
                    }
                    result.push_back(result_row);
                    result_row = {};
                }
                return Matrix(result);
            } else {
                throw std::invalid_argument("Матрицы не совпадают по размеру");
            }
        } else {
            throw std::invalid_argument("Матрицы не совпадают по размеру");;
        }
    }
//_____________________________________________________________________________________________________________________
    //Эта функция будет отвечать за разность двух матриц.
    Matrix dif_mtx(const Matrix& other) const {
        // собираем данные для того чтобы не вызывать функцию повторно
        int row1 = data.size();
        int row2 = other.data.size();
        int coloms1 = data[0].size();
        int coloms2 = other.data[0].size();
        // Также создаём два вектора, один из которых мы будем использовать как строки другого,
        // а второй будет результатом выполнения функции
        std::vector<std::vector<double>> result;
        std::vector<double> result_row;


        // Переходим к проверке условий и реализации
        if (row1 == row2){
            if (coloms1 == coloms2){
                for (int i = 0; i < row1; i++){
                    for (int j = 0; j < coloms1; j++){
                        double elem = data[i][j] - other.data[i][j];
                        result_row.push_back(elem);
                    }
                    result.push_back(result_row);
                    result_row = {};
                }
                return Matrix(result);
            } else {
                throw std::invalid_argument("Матрицы не совпадают по размеру");
            }
        } else {
            throw std::invalid_argument("Матрицы не совпадают по размеру");;
        }
    }
//_____________________________________________________________________________________________________________________
    // реализация унарного минуса (знак всех элементов умножается на противоположный)
    void unary_minus(){
        for (int i = 0; i < data.size(); i++){
            for (int j = 0; j < data[0].size(); j++){
                data[i][j] *= -1;
            }
        }
    }
//_____________________________________________________________________________________________________________________
    // Умножение матрицы на число
    Matrix product_num(double num){
        std::vector<std::vector<double>> result (data.size(), std::vector<double>(data[0].size(), 1));
        for(int i = 0; i < data.size(); i++){
            for (int j = 0; j < data[0].size(); j++){
                double tmp = data[i][j] * num;
                result[i][j] *= tmp;
            }
        }
        return Matrix(result);
    }
//_____________________________________________________________________________________________________________________
    // Функция перемножения матрицы
    // Если первая матрица будет размером m * n, а вторая размером n * p, то конечная матрица
    // будет размером m * p
    Matrix product(const Matrix& other) const{
        int row1 = data.size();
        int coloms2 = other.data[0].size();
        int n = other.data.size();
        // Также создаём переменные с данными
        // Необычным стало то, что тут пришлось задать вектор необходимых размеров и заполнить его нулями
        // возможно другое решение, но именно это решение даёт возможность не менять алгоритм
        std::vector<std::vector<double>> result(row1, std::vector<double>(coloms2, 0));

        if (data[0].size() == other.data.size()) {
            for (int i = 0; i < row1; i++) {
                for (int j = 0; j < coloms2; j++) {
                    for (int k = 0; k < n; k++) {
                        result[i][j] += data[i][k] * other.data[k][j];
                    }
                }
            }
            return Matrix(result);
        } else {
            throw std::invalid_argument("Матрицы не подходя по размеру(кол-во строк во второй не равно "
                                        "кол-ву столбцов в первой)");
        }
    }
//_____________________________________________________________________________________________________________________
    // На этом этапе я сталкнулся с интересной делеммой: чтобы удобно дилить матрицу на матрицу, мне удобно искать
    // обратную матрицу и умножать, но чтобы искать обратную матрицу, нужно реализовать деление матрицы на константу,
    // а я планировал сделать это только после реализации деления матрицы на матрицу ))
    //Функция деления матрицы на константу (число)
    Matrix division_num(double num){
        std::vector<std::vector<double>> result(data.size(), std::vector<double>(data[0].size(), 1));
        for(int i = 0; i < data.size(); i++){
            for (int j = 0; j < data[0].size(); j++){
                double tmp = data[i][j] / num;
                result[i][j] = tmp;
            }
        }
        return Matrix(result);
    }
//_____________________________________________________________________________________________________________________
    // Теперь я могу реализовать ФУНКЦИЮ ПОЛУЧЕНИЯ ОБРАТНОЙ МАТРИЦЫ
    Matrix inverse() const{
    // Для начала нужно проверить, не нулевая ли матрица и существует ли она
    double det = determinant();
    if (det != 0){
        int n = data.size();
        std::vector<std::vector<double>> result(n, std::vector<double>((n), 1));
        // Понимаю что уже не важно data.size() или data[0].size(), так что легче ввести это значение в переменную,
        // но я не делаю так весьде потому что потом мне будет трудно связать переменные с их значениями
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double sign = (i + j) % 2 == 0 ? 1.0 : -1.0;
                result[j][i] = sign * cofactor(i, j);
                }
            }
        // Такой странный способ сделан для того, чтобы связать функцию деления на константу с функцией получения
        // обратной матрицы (легче было просто скопировать алгоритм)
        Matrix final_data;
        final_data = Matrix(result).division_num(det);
        return final_data;
        }
    else { // Вырожденной матрицой называют матрицу с определителем == 0
        throw std::invalid_argument("Матрица вырожденная или у неё нет определителя");;
        }

    }

//_____________________________________________________________________________________________________________________
    // Наконец я дошел до этого
    // Функция деления матрицы на матрицу!
    Matrix divicion(const Matrix& other) const{
        return product(other.inverse());
    }
//_____________________________________________________________________________________________________________________
    //Функция вычисления ранга матрицы
    int rank() const {
        int n = data.size();    // Количество строк
        int m = data[0].size(); // Количество столбцов
        std::vector<std::vector<double>> a = data; // Копия данных матрицы для преобразований

        int rank = 0;
        for (int col = 0; col < m && rank < n; ++col) {
            // Ищем ненулевой элемент в текущем столбце
            int swap_row = rank;
            while (swap_row < n && fabs(a[swap_row][col]) < 1e-10) {
                swap_row++;
            }
            if (swap_row == n) {
                continue; // Все элементы в текущем столбце равны 0
            }

            // Переставляем строки, если нашли ненулевой элемент
            if (swap_row != rank) {
                std::swap(a[swap_row], a[rank]);
            }

            // Делаем элементы текущего столбца ниже a[rank][col] равными 0
            for (int row = rank + 1; row < n; ++row) {
                double coef = a[row][col] / a[rank][col];
                for (int j = col; j < m; ++j) {
                    a[row][j] -= a[rank][j] * coef;
                }
            }

            // Увеличиваем ранг, так как нашли линейно независимую строку
            ++rank;
        }

        return rank;
    }
//_____________________________________________________________________________________________________________________
    Matrix transpose() const {

        std::vector<std::vector<double>> transposedData(data[0].size(), std::vector<double>(data.size()));

        for (int i = 0; i < data.size(); i++) {
            for (int j = 0; j < data[0].size(); j++) {
                transposedData[j][i] = data[i][j];
            }
        }

        return Matrix(transposedData);
    }
//_____________________________________________________________________________________________________________________
    // LU-разложение
    std::pair<Matrix, Matrix> luDecomposition() const {
        int n = data.size();
        std::vector<std::vector<double>> l(n, std::vector<double>(n, 0.0));
        std::vector<std::vector<double>> u(data);

        for (int i = 0; i < n; i++) {
            l[i][i] = 1.0; // Диагональные элементы L равны 1
            for (int j = i + 1; j < n; j++) {
                if (u[i][i] == 0) throw std::runtime_error("Элемент матрицы равен 0, деление невозможно");

                double factor = u[j][i] / u[i][i];
                l[j][i] = factor;
                for (int k = i; k < n; k++) {
                    u[j][k] -= factor * u[i][k];
                }
            }
        }

        return {Matrix(l), Matrix(u)};
    }
//_____________________________________________________________________________________________________________________
    //QR-разложение (с использованием метода Грама-Шмидта)
    std::pair<Matrix, Matrix> qrDecomposition() const {
        std::vector<std::vector<double>> Q(data.size(), std::vector<double>(data[0].size(), 0));
        std::vector<std::vector<double>> R(data[0].size(), std::vector<double>(data[0].size(), 0));

        // Векторы для хранения ортогональных векторов
        std::vector<std::vector<double>> u(data[0].size(), std::vector<double>(data.size(), 0));

        for (int i = 0; i < data[0].size(); i++) {
            // Получаем i-й столбец матрицы data
            for (int j = 0; j < data.size(); j++) {
                u[i][j] = data[j][i];
            }

            // Ортогонализация с помощью метода Грама-Шмидта
            for (int j = 0; j < i; j++) {
                double dot = 0;
                for (int k = 0; k < data.size(); k++) {
                    dot += u[i][k] * Q[k][j];
                }
                R[j][i] = dot;
                for (int k = 0; k < data.size(); k++) {
                    u[i][k] -= dot * Q[k][j];
                }
            }

            // Нормализация вектора u[i] и запись в Q
            double norm = 0;
            for (int j = 0; j < data.size(); j++) {
                norm += u[i][j] * u[i][j];
            }
            norm = sqrt(norm);
            R[i][i] = norm;
            for (int j = 0; j < data.size(); j++) {
                Q[j][i] = u[i][j] / norm;
            }
        }

        return {Matrix(Q), Matrix(R)};
    }
//_____________________________________________________________________________________________________________________
};


int main() {
    Matrix m5({
                      {4, 3, 2},
                      {3, 5, 1},
                      {2, 1, 3}
              });
    auto [L, U] = m5.luDecomposition();
    L.print();
    U.print();

    Matrix m6({
                      {12, -51, 4},
                      {6, 167, -68},
                      {-4, 24, -41}
    });
    auto [Q, R] = m6.qrDecomposition();
    std::cout << "Ортогональная матрица Q" << std::endl;
    Q.print();
    std::cout << "Верхняя треугольная матрица R" << std::endl;
    R.print();

    Matrix m1(3, 3, 3);
    std::vector<std::vector<double>> vec = {
            {1, 2, 4, 7, 11},
            {2, 3, 5, 8, 12},
            {3, 4, 6, 9, 13}
    };
    Matrix m2(vec);
    Matrix m3(3, 3, 1);

    m1.print();
    m1.unary_minus();
    std::cout << m1.rank() << std::endl;
    m1.print();
    m2.print();
    Matrix m2_tr = m2.transpose();
    m2_tr.print();

    std::vector<std::vector<double>> cat_face = {
            {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0},
            {0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0},
            {0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0},
            {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
            {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
            {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
            {0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0},
            {0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0}
    };
    Matrix m_cat(cat_face);
    m_cat.print_cat();

    std::cout << m1.get_coloms() << " ";
    std::cout << m2.find_elem(2, 4) << std::endl;
    std::cout<< m2.find_elem(0, 0) << std::endl;

    Matrix sum_1 = m1.sum_mtx(m3);
    Matrix dif_1 = m1.dif_mtx(m3);
    sum_1.print();
    dif_1.print();

    Matrix pr_1 = m1.product(m3);
    pr_1.print();
    double h = 7;
    Matrix pr_c_1 = m1.product_num(h);
    pr_c_1.print();

    return 0;
}
// Вывод: мне очень сильно надоело соблюдать табуляцию и искать решение к каждой проблеме, которая выскакивает.
// Изначально я подумал, что задание очень простое, рисовал рисунки, сейчас понимаю что дойти до тестов будет
// для меня личной победой. Я снова хочу увидеть котика, которого рисует мне матрица в комментариях.

//TEST_CASE("Testing_Matrix"){
//    Matrix m1(3, 3, 3);
//    std::vector<std::vector<double>> vec = {
//            {1, 2, 4, 7, 11},
//            {2, 3, 5, 8, 12},
//            {3, 4, 6, 9, 13}
//    };
//    Matrix m2(vec);
//    Matrix m3 = m1.product(m2);
//
//    CHECK(m1.determinant() == 0);
//    CHECK(m2.get_row() == 3);
//    CHECK(m2.get_coloms() == 5);
//    CHECK(m2.find_elem(2, 1) == 4);
//
//    m2.set_elem(2, 1, 10);
//
//    CHECK(m2.find_elem(2, 1) == 10);
//
//
//    Matrix m5({
//                  {4, 3, 2},
//                  {3, 5, 1},
//                  {2, 1, 3}
//          });
//    Matrix l1({
//        {1, 0, 0},
//        {0.75, 1, 0},
//        {0.5, -0.181818, 1}
//    });
//    Matrix u1({
//                      {4, 3, 2},
//                      {0, 2.75, -0.5},
//                      {0, 1.38778e-17, 1.90909}
//    });
//    [auto L,auto U] = m5.luDecomposition();
//    CHECK([auto L,auto U] == [l1, u1]);
//}