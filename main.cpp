#include <iostream>
#include <vector>
#include <stdexcept>
#include <cmath>

template <typename T>

class matrix{
private:
    size_t rows;
    size_t cols;
    std::vector<std::vector<T>> data;
public:
    matrix() : rows(0), cols(0) {}

    matrix(size_t r, size_t c, T InitValue = T()) : rows(r), cols(c), data(r, std::vector<T>(c, InitValue)) {}

    matrix(const matrix& other) : rows(other.rows), cols(other.cols), data(other.data) {}

    ~matrix() {}

    matrix inverse() const{
        if (rows != cols) throw std::invalid_argument("For inverse rows matrix must = cols matrix");

        size_t n = rows;
        matrix doublematrix(n, 2*n);
        for(size_t i = 0; i < n; ++i){
            for(size_t j = 0; j < n; ++j){
                doublematrix.data[i][j] = data[i][j];
            }
            doublematrix.data[i][n + i] = 1;
        }

        for(size_t i = 0; i < n; ++i){
            T pivot = doublematrix.data[i][i];
            if(std::abs(pivot) = 1e-9) throw std::runtime_error("Matrix is singular and cannot be inverted.");

            for(size_t j = 0; j < 2*n; ++j){
                doublematrix.data[i][j] /= pivot;
            }
            for(size_t k = 0; k < n; ++k){
                if(k != i){
                    T factor = doublematrix.data[k][i];
                    for(size_t j = 0; j < 2*n; ++j){
                        doublematrix.data[k][j] -= factor * doublematrix.data[i][j];
                    }
                }
            }
        }

        matrix invers(n, n);
        for(size_t i = 0; i < n; ++i){
            for(size_t j = 0; j < n; ++j){
                invers.data[i][j] = doublematrix.data[i][j + n];
            }
        }
        return invers;

    }

    // ======================================================================================================================
    matrix& operator = (const matrix& other){
        if(this != &other){
            rows = other.rows;
            cols = other.cols;
            data = other.data;
        }
        return *this;
    }

    bool operator == (const matrix& other) const{ // показываем, что мы ничего не делаем с объектом
        return this->data == other.data;
    }

    bool operator != (const matrix& other) const{
        return !(*this == other);
    }
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    matrix operator + (const matrix& other) const{
        if(rows != other.rows || cols != other.cols) throw std::invalid_argument("\nfor operator +, the dimensions of the matrices must match\n");
        matrix result(rows, cols);
        for(size_t i = 0; i < rows; ++i){
            for(size_t j = 0; j < cols; ++j){
                result.data[i][j] = data[i][j] + other.data[i][j];
            }
        }
        return result;
    }

    matrix operator + (T scalar) const{
        matrix result(rows, cols);
        for(size_t i = 0; i < rows; ++i){
            for(size_t j = 0; j < cols; ++j){
                result.data[i][j] = data[i][j] + scalar;
            }
        }
        return result;
    }

    friend matrix operator + (T scalar, const matrix& other){
        return other + scalar;
    }


    matrix& operator += (const matrix& other){
        if(rows != other.rows || cols != other.cols) throw std::invalid_argument("\nfor operator +, the dimensions of the matrices must match\n");
        for(size_t i = 0; i < rows; ++i){
            for(size_t j = 0; j < cols; ++j){
                data[i][j] += other.data[i][j];
            }
        }
        return *this;
    }

    matrix& operator ++ (){
        for(auto& row : data)
            for(auto& item : row)
                ++item;
        return this;
    }

    matrix operator ++ (int){
        matrix temp(*this); // чтобы вернуть после инкрементации(увеличении элемента на единицу)
        ++(*this);
        return temp;
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------
    matrix operator - (const matrix& other) const{
        if(rows != other.rows || cols != other.cols) throw std::invalid_argument("\nfor operator -, the dimensions of the matrices must match\n");
        matrix result(rows, cols);
        for(size_t i = 0; i < rows; ++i){
            for(size_t j = 0; j < cols; ++j){
                result.data[i][j] = data[i][j] - other.data[i][j];
            }
        }
        return result;
    }

    matrix operator - (T scalar) const{
        matrix result(rows, cols);
        for(size_t i = 0; i < rows; ++i){
            for(size_t j = 0; j < cols; ++j){
                result.data[i][j] = data[i][j] - scalar;
            }
        }
        return result;
    }

    friend matrix operator - (T scalar, const matrix& other){
        matrix result(other.rows, other.cols);
        for(size_t i = 0; i < other.rows; ++i){
            for(size_t j = 0; j < other.cols; ++j){
                result.data[i][j] = scalar - other.data[i][j];
            }
        }
        return result;
    }


    matrix& operator -= (const matrix& other){
        if(rows != other.rows || cols != other.cols) throw std::invalid_argument("\nfor operator -, the dimensions of the matrices must match\n");
        for(size_t i = 0; i < rows; ++i){
            for(size_t j = 0; j < cols; ++j){
                data[i][j] -= other.data[i][j];
            }
        }
        return *this;
    }

    matrix& operator -- (){
        for(auto& row : data)
            for(auto& item : row)
                --item;
        return this;
    }

    matrix operator -- (int){
        matrix temp(*this); // чтобы вернуть после инкрементации(увеличении элемента на единицу)
        --(*this);
        return temp;
    }
    //********************************************************************************************************************** */
    matrix operator * (const matrix& other) const{
        if(cols != other.rows) throw std::invalid_argument("\nfor operator *, cols matrix 1 must = rows matrix 2\n");
        matrix result(rows, other.cols);

        for(size_t i = 0; i < rows; ++i){
            for(size_t j = 0; j < other.cols; ++j){
                result.data[i][j] = T();
                for(size_t k = 0; k < cols; ++k){
                    result.data[i][j] += data[i][k] * other.data[k][j];
                }
            }
        }
        return result;
    }

    matrix operator * (T scalar) const{
        matrix result(rows, cols);
        for(size_t i = 0; i < rows; ++i){
            for(size_t j = 0; j < cols; ++j){
                result.data[i][j] = data[i][j] * scalar;
            }
        }
        return result;
    }

    friend matrix operator * (T scalar, const matrix& other){
        return other * scalar;
    }


    matrix& operator *= (const matrix& other){
        if(cols != other.rows) throw std::invalid_argument("\nfor operator *, cols matrix 1 must = rows matrix 2\n");
        matrix result(rows, other.cols);

        for(size_t i = 0; i < rows; ++i){
            for(size_t j = 0; j < other.cols; ++j){
                result.data[i][j] = T();
                for(size_t k = 0; k < cols; ++k){
                    result.data[i][j] += data[i][k] * other.data[k][j];
                }
            }
        }
        *this = result;
        return *this;
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    matrix operator / (const matrix& other) const{
        return *this * other.inverse();
    }

    matrix operator / (T scalar) const{
        matrix result(rows, cols);
        for(size_t i = 0; i < rows; ++i){
            for(size_t j = 0; j < cols; ++j){
                result.data[i][j] = data[i][j] / scalar;
            }
        }
        return result;
    }

    friend matrix operator / (T scalar, const matrix& other){
        matrix result(other.rows, other.cols);
        for(size_t i = 0; i < other.rows; ++i){
            for(size_t j = 0; j < other.cols; ++j){
                result.data[i][j] = scalar / other.data[i][j];
            }
        }
        return result;
    }


    matrix& operator /= (const matrix& other){
        if(cols != other.rows) throw std::invalid_argument("\nfor operator *, cols matrix 1 must = rows matrix 2\n");
        matrix result(rows, other.cols);
        //tytyytytyt
        return result;
    }
    //()((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((()))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))

    T& operator()(size_t i, size_t j){
        if(i >= rows) throw std::out_of_range("matrix index out of rows");
        if(j >= cols) throw std::out_of_range("matrix index out of cols");
        return data[i][j];
    }

    const T& operator()(size_t i, size_t j) const{
        if(i >= rows) throw std::out_of_range("matrix index out of rows");
        if(j >= cols) throw std::out_of_range("matrix index out of cols");
        return data[i][j];
    }

    // prepre

    explicit operator size_t() const{
        return rows * cols;
    }

    // func

    void print() const{
        for(const auto& row : data){
            for(const auto& item : row){
                std::cout << item << " ";
            }
            std::cout << std::endl;
        }
    }
};

int main(){
    matrix goida(3, 3, 3);
    matrix baza(3, 3, 52);
    printf("Matrix 1\n");
    baza.print();
    printf("\n");
    printf("Matrix 2\n");
    goida.print();

    // операторы присваивания
    printf("\nOPERATORS != = ==\n");
    printf("Matrix 1 == Matrix 2 : %d\n", (int) (baza==goida));
    goida = baza;
    printf("Matrix 1 = Matrix 2\n");
    printf("Matrix 1 == Matrix 2 : %d\n", (int) (baza!=goida));

    // операторы + - * /
    printf("\n+ OPERATORS +\n");
    printf("Matrix 1\n");
    baza.print();
    printf("\n");
    printf("Matrix 2\n");
    goida.print();

    matrix temp(3, 3, 0);
    temp = baza + goida;
    printf("\nMatrix 1 + Matrix 2 : \n");
    temp.print();
    temp = goida + 5;
    printf("\nMatrix 1 + 5 : \n");
    temp.print();
    temp = 10 + goida;
    printf("\n10 + Matrix 1 : \n");
    temp.print();

    //----------------------------

    printf("\n- OPERATORS -\n");
    printf("Matrix 1\n");
    baza.print();
    printf("\n");
    printf("Matrix 2\n");
    goida.print();

    temp = baza - goida;
    printf("\nMatrix 1 - Matrix 2 : \n");
    temp.print();
    temp = goida - 5;
    printf("\nMatrix 1 + 5 : \n");
    temp.print();
    temp = 100 - goida;
    printf("\n100 - Matrix 1 : \n");
    temp.print();

    //************************************* */

    printf("\n* OPERATORS *\n");
    printf("Matrix 1\n");
    baza.print();
    printf("\n");
    printf("Matrix 2\n");
    goida.print();

    temp = baza * goida;
    printf("\nMatrix 1 * Matrix 2 : \n");
    temp.print();
    temp = goida * 5;
    printf("\nMatrix 1 * 5 : \n");
    temp.print();
    temp = 100 / goida;
    printf("\n100 * Matrix 1 : \n");
    temp.print();

    ///////////////////////////////////////////////////////////
    printf("\n/ OPERATORS /\n");
    printf("Matrix 1\n");
    baza.print();
    printf("\n");
    printf("Matrix 2\n");
    goida.print();

    //temp = baza / goida;
    //printf("\nMatrix 1 / Matrix 2 : \n");
    //temp.print();
    temp = goida / 5;
    printf("\nMatrix 1 / 5 : \n");
    temp.print();
    temp = 104 / goida;
    printf("\n104 / Matrix 1 : \n");
    temp.print();

    //Арифметика с накоплением (+=,-=);
    printf("\n\n ============================= \n\n");
    //+=+=+=+=+=
    printf("\n+= OPERATORS +=\n");
    printf("Matrix 1\n");
    baza.print();
    printf("\n");
    printf("Matrix 2\n");
    goida.print();

    goida += baza;
    printf("\nMatrix 2 += Matrix 1\n");
    goida.print();

    //------------------------------

    printf("\n-= OPERATORS -=\n");
    printf("Matrix 1\n");
    baza.print();
    printf("\n");
    printf("Matrix 2\n");
    goida.print();

    goida -= baza;
    printf("\nMatrix 2 -= Matrix 1\n");
    goida.print();

    //********************************************** */
    printf("\n*= OPERATORS *=\n");
    printf("Matrix 1\n");
    baza.print();
    printf("\n");
    printf("Matrix 2\n");
    goida.print();

    goida *= baza;
    printf("\nMatrix 2 *= Matrix 1\n");
    goida.print();
    //Унарные (++,--) в префиксной и постфиксной форме;
    //Логические (<, >, ==, != );
    //Операторы взятия элемента ( [] или () ) по номеру или ключу;
    //Операторы преобразования типа к любому базовому.
}