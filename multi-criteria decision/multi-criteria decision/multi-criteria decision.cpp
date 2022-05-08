#include <iostream>
#include <vector>
#include <string>

int readPositiveInt(std::string input_message, int upper_bound = INT_MAX) {
    int result;
    for (;;) {
        std::cout << input_message << std::flush;
        if ((std::cin >> result).good() && (1 <= result) && (result <= upper_bound))
            return result;
        if (std::cin.fail()) {
            std::cin.clear();
            std::cout << "Неверный ввод.\n";
        } else {
            std::cout << "Число должно быть положительным и меньше " << upper_bound << "\n";
        }
        std::cin.ignore(std::numeric_limits<std::size_t>::max(), '\n');
    }
}

double readPositiveDouble(double upper_bound = 1) {
    double result;
    for (;;) {
        if ((std::cin >> result).good() && (0 <= result) && (result <= upper_bound + 0.001))
            return result;
        if (std::cin.fail()) {
            std::cin.clear();
            std::cout << "Некорректный ввод.\n";
        } else {
            std::cout.precision(3);
            std::cout << "Число должно быть положительным и меньше " << upper_bound <<"\n";
        }
        std::cin.ignore(std::numeric_limits<std::size_t>::max(), '\n');
    }
}

double sum(std::vector<double>* values) {
    double sum = 0;
    for (auto val : *values) {
        sum += val;
    }
    return sum;
}

void inputWeights(std::string input, std::vector<double>* weights) {
    do {
        std::cout << input;
        for (size_t i = 0; i < weights->size(); ++i) {
            (*weights)[i] = readPositiveDouble();
        }
    } while (abs(sum(weights) - 1) > 0.0001);
}

double weightedAverage(std::vector<double>* values, std::vector<double>* weights, int size) {
    double result = 1;
    for (int i = 0; i < size; ++i) {
        result *= pow((*values)[i],(*weights)[i]);
    }
    return result;
}

std::vector<double> getMinChoice(std::vector<double>* weights, std::vector<std::vector<double>>* values,
                                 int size, int alternative_c) {
    std::vector<double> min_values(alternative_c);
    double temp;
    for (size_t i = 0; i < alternative_c; ++i) {
        min_values[i] = pow((*values)[i][0], (*weights)[0]);
        for (int j = 1; j < size; ++j) {
            temp = pow((*values)[i][j], (*weights)[j]);
            if (temp < min_values[i]) {
                min_values[i] = temp;
            }
        }
    }
    return min_values;
}

std::pair<double, double> intuitWeightedAverage(std::vector<double>* mu_values, std::vector<double>* nu_values,
                                                std::vector<double>* weights, int size) {
    double mu = 1, nu = 1;
    for (int i = 0; i < size; ++i) {
        mu *= pow((1 - (*mu_values)[i]), (*weights)[i]);
        nu *= pow((*nu_values)[i], (*weights)[i]);
    }
    return std::make_pair(mu, nu);
}

std::vector<double> inputCriteriaWeights(int criteria_count, int expert_count, std::vector<double>* expert_weights) {
    std::vector<std::vector<double>> expert_criteria_weights(criteria_count, std::vector<double>(expert_count));
    for (int i = 0; i < criteria_count; ++i) {
        for (int j = 0; j < expert_count; ++j) {
            std::cout << "Введите вес, которые определил " + std::to_string(j + 1)
                      << " эксперт для критерия " + std::to_string(i + 1) + " ";
            expert_criteria_weights[i][j] = readPositiveDouble();
        }
    }

    std::vector<double> criteria_weights(criteria_count);
    for (int i = 0; i < criteria_count; ++i) {
        criteria_weights[i] = weightedAverage(&expert_criteria_weights[i], expert_weights, criteria_count);
    }
    return criteria_weights;
}

int main() {
    setlocale(LC_ALL, "rus");
    int expert_count = readPositiveInt("Введите число экспертов: ");
    int alternative_count = readPositiveInt("Введите число альтернатив: ");
    int criteria_count = readPositiveInt("Введите число критериев: ");
    std::vector<double> expert_weights(expert_count);
    inputWeights("Введите веса экспертов:\n", &expert_weights);
    std::vector<double> criteria_weights = inputCriteriaWeights(criteria_count, expert_count, &expert_weights);

    std::vector<double> choices;
    int choice = readPositiveInt("КНМ(1) или ИНМ(2)? ", 2);
    if (choice == 1) {
        int count;
        double sum = 0;
        std::vector<std::vector<double>> s_h_elem(criteria_count, std::vector<double>(alternative_count));
        for (int i = 0; i < criteria_count; ++i) {
            for (int j = 0; j < alternative_count; ++j) {
                sum = 0;
                count = readPositiveInt("Количество значений функции принадлежности для " +
                                        std::to_string(j + 1) + " альтернативы " +
                                        std::to_string(i + 1) + " критерию ", expert_count * 2);
                for (int k = 0; k < count; ++k) {
                    std::cout << "Значение " << k + 1 << " = ";
                    sum += readPositiveDouble();
                }
                s_h_elem[i][j] = sum / count;
            }
        }
        choices = getMinChoice(&criteria_weights, &s_h_elem, criteria_count, alternative_count);

    } else {
        std::pair<double, double> aggreg;
        std::vector<std::vector<double>> w_func(criteria_count, std::vector<double>(alternative_count));
        std::vector<double> expert_mu(expert_count);
        std::vector<double> expert_nu(expert_count);
        for (int i = 0; i < criteria_count; ++i) {
            for (int j = 0; j < alternative_count; ++j) {
                for (int k = 0; k < expert_count; ++k) {
                    std::cout << "Значение функции принадлежности " + std::to_string(i + 1) + " альтернативы "
                              << std::to_string(j + 1) + " критерию " + std::to_string(k + 1) + " эксперта ";
                    expert_mu[k] = readPositiveDouble();
                    std::cout << "Значение функции непринадлежности " + std::to_string(i + 1) + " альтернативы "
                              << std::to_string(j + 1) + " критерию " + std::to_string(k + 1) + " эксперта ";
                    expert_nu[k] = readPositiveDouble(1 -  expert_mu[k]);
                }
                aggreg = intuitWeightedAverage(&expert_mu, &expert_nu, &expert_weights, expert_count);
                w_func[i][j] = 0.5 * aggreg.first + 1.5 * (1 - aggreg.second) - 1;
            }
        }
        choices = getMinChoice(&criteria_weights, &w_func, criteria_count, alternative_count);
    
    }
    std::cout.precision(3);
    double max = choices[0];
    int index_max = 0;
    for (int i = 0; i < criteria_count; ++i) {
        std::cout << "Альтернатива " << i + 1 << ": " << choices[i] << '\n';
        if (choices[i] > max) {
            index_max = i;
            max = choices[i];
        }
    }
    std::cout << "Наиболее предпочтительная альтернатива - " << index_max + 1;
    return 0;
}
