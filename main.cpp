#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <cmath>
#include <fstream>
#include <set>
#include <tuple>

const int Nd = 9;

class Candidate{
public:
    std::vector<std::vector<int>> values;
    double fitness;

    Candidate() : values(Nd, std::vector<int>(Nd, 0)), fitness(0.0) {}

    void update_fitness();
    bool mutate(double mutation_rate, const Candidate &given);
    bool is_row_duplicate(int row, int value) const;
    bool is_column_duplicate(int column, int value) const;
    bool is_block_duplicate(int row, int column, int value) const;
};

class Population{
public:
    std::vector<Candidate> candidates;

    void seed(int Nc, const Candidate &given);
    void update_fitness();
    void sort();
    static bool sort_fitness(const Candidate &x, const Candidate &y);
};

class Given : public Candidate{
public:
    Given(const std::vector<std::vector<int>> &values){
        this->values = values;
    }
};

class Tournament{
public:
    static Candidate compete(const std::vector<Candidate> &candidates);
};

class CycleCrossover{
public:
    static std::pair<Candidate, Candidate> crossover(const Candidate &parent1, const Candidate &parent2, double crossover_rate);
    static std::pair<std::vector<int>, std::vector<int>> crossover_rows(const std::vector<int> &row1, const std::vector<int> &row2);
    static int find_unused(const std::vector<int> &parent_row, const std::vector<int> &remaining);
    static int find_value(const std::vector<int> &parent_row, int value);
};

class Sudoku{
public:
    Candidate given;
    Population population;

    void load(const std::string &path);
    void save(const std::string &path, const Candidate &solution);
    Candidate solve();
};

// Function implementations

void Candidate::update_fitness(){
    std::vector<int> row_count(Nd, 0), column_count(Nd, 0), block_count(Nd, 0);
    double row_sum = 0.0, column_sum = 0.0, block_sum = 0.0;

    for (int i = 0; i < Nd; ++i){
        for (int j = 0; j < Nd; ++j){
            row_count[values[i][j] - 1]++;
        }
        row_sum += (1.0 / std::set<int>(row_count.begin(), row_count.end()).size()) / Nd;
        std::fill(row_count.begin(), row_count.end(), 0);
    }

    for (int i = 0; i < Nd; ++i){
        for (int j = 0; j < Nd; ++j){
            column_count[values[j][i] - 1]++;
        }
        column_sum += (1.0 / std::set<int>(column_count.begin(), column_count.end()).size()) / Nd;
        std::fill(column_count.begin(), column_count.end(), 0);
    }

    for (int i = 0; i < Nd; i += 3){
        for (int j = 0; j < Nd; j += 3){
            for (int k = 0; k < 3; ++k){
                for (int l = 0; l < 3; ++l){
                    block_count[values[i + k][j + l] - 1]++;
                }
            }
            block_sum += (1.0 / std::set<int>(block_count.begin(), block_count.end()).size()) / Nd;
            std::fill(block_count.begin(), block_count.end(), 0);
        }
    }

    if (static_cast<int>(row_sum) == 1 && static_cast<int>(column_sum) == 1 && static_cast<int>(block_sum) == 1){
        fitness = 1.0;
    }else{
        fitness = column_sum * block_sum;
    }
}

bool Candidate::mutate(double mutation_rate, const Candidate &given){
    static std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<double> dist(0.0, 1.1);
    double r = dist(rng);

    while (r > 1.0){
        r = dist(rng);
    }

    if (r < mutation_rate){
        bool success = false;
        while (!success){
            int row1 = std::uniform_int_distribution<int>(0, 8)(rng);
            int from_column = std::uniform_int_distribution<int>(0, 8)(rng);
            int to_column = std::uniform_int_distribution<int>(0, 8)(rng);

            while (from_column == to_column){
                from_column = std::uniform_int_distribution<int>(0, 8)(rng);
                to_column = std::uniform_int_distribution<int>(0, 8)(rng);
            }

            if (given.values[row1][from_column] == 0 && given.values[row1][to_column] == 0){
                if (!given.is_column_duplicate(to_column, values[row1][from_column]) &&
                    !given.is_column_duplicate(from_column, values[row1][to_column]) &&
                    !given.is_block_duplicate(row1, to_column, values[row1][from_column]) &&
                    !given.is_block_duplicate(row1, from_column, values[row1][to_column])){
                    std::swap(values[row1][from_column], values[row1][to_column]);
                    success = true;
                }
            }
        }
        return true;
    }
    return false;
}

void Population::seed(int Nc, const Candidate &given){
    candidates.clear();
    Candidate helper;
    helper.values.resize(Nd, std::vector<int>(Nd, 0));

    std::vector<std::vector<std::vector<int>>> possible_values(Nd, std::vector<std::vector<int>>(Nd));

    for (int row = 0; row < Nd; ++row){
        for (int column = 0; column < Nd; ++column){
            for (int value = 1; value <= 9; ++value){
                if (given.values[row][column] == 0 &&
                    !given.is_column_duplicate(column, value) &&
                    !given.is_block_duplicate(row, column, value) &&
                    !given.is_row_duplicate(row, value)){
                    possible_values[row][column].push_back(value);
                } else if (given.values[row][column] != 0){
                    possible_values[row][column].push_back(given.values[row][column]);
                    break;
                }
            }
        }
    }

    std::random_device rd;
    std::mt19937 gen(rd());

    for (int p = 0; p < Nc; ++p){
        Candidate g;
        for (int i = 0; i < Nd; ++i){
            std::vector<int> row(Nd);
            for (int j = 0; j < Nd; ++j){
                if (given.values[i][j] != 0){
                    row[j] = given.values[i][j];
                } else {
                    std::uniform_int_distribution<> dis(0, possible_values[i][j].size() - 1);
                    row[j] = possible_values[i][j][dis(gen)];
                }
            }

            while (std::set<int>(row.begin(), row.end()).size() != Nd){
                for (int j = 0; j < Nd; ++j){
                    if (given.values[i][j] == 0){
                        std::uniform_int_distribution<> dis(0, possible_values[i][j].size() - 1);
                        row[j] = possible_values[i][j][dis(gen)];
                    }
                }
            }
            g.values[i] = row;
        }
        candidates.push_back(g);
    }
    update_fitness();
    std::cout << "Seeding complete." << std::endl;
}

void Population::update_fitness(){
    for (auto &candidate : candidates){
        candidate.update_fitness();
    }
}

void Population::sort(){
    std::sort(candidates.begin(), candidates.end(), sort_fitness);
}

bool Population::sort_fitness(const Candidate &x, const Candidate &y){
    return x.fitness > y.fitness;
}

bool Candidate::is_row_duplicate(int row, int value) const{
    for (int column = 0; column < Nd; ++column){
        if (values[row][column] == value){
            return true;
        }
    }
    return false;
}

bool Candidate::is_column_duplicate(int column, int value) const {
    for (int row = 0; row < Nd; ++row){
        if (values[row][column] == value){
            return true;
        }
    }
    return false;
}

bool Candidate::is_block_duplicate(int row, int column, int value) const{
    int i = 3 * (row / 3);
    int j = 3 * (column / 3);

    for (int r = 0; r < 3; ++r){
        for (int c = 0; c < 3; ++c){
            if (values[i + r][j + c] == value){
                return true;
            }
        }
    }
    return false;
}

Candidate Tournament::compete(const std::vector<Candidate> &candidates){
    static std::mt19937 rng(std::random_device{}());
    int index1 = std::uniform_int_distribution<int>(0, candidates.size() - 1)(rng);
    int index2 = std::uniform_int_distribution<int>(0, candidates.size() - 1)(rng);

    const Candidate &c1 = candidates[index1];
    const Candidate &c2 = candidates[index2];

    const Candidate &fittest = (c1.fitness > c2.fitness) ? c1 : c2;
    const Candidate &weakest = (c1.fitness > c2.fitness) ? c2 : c1;

    double selection_rate = 0.85;
    std::uniform_real_distribution<double> dist(0.0, 1.1);
    double r = dist(rng);

    while (r > 1.0){
        r = dist(rng);
    }

    return (r < selection_rate) ? fittest : weakest;
}

std::pair<Candidate, Candidate> CycleCrossover::crossover(const Candidate &parent1, const Candidate &parent2, double crossover_rate){
    Candidate child1, child2;

    child1.values = parent1.values;
    child2.values = parent2.values;

    static std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<double> dist(0.0, 1.1);
    double r = dist(rng);

    while (r > 1.0){
        r = dist(rng);
    }

    if (r < crossover_rate){
        std::uniform_int_distribution<int> dist(0, 8);
        int crossover_point1 = dist(rng);
        int crossover_point2 = dist(rng);

        while (crossover_point1 == crossover_point2){
            crossover_point2 = dist(rng);
        }

        if (crossover_point1 > crossover_point2){
            std::swap(crossover_point1, crossover_point2);
        }

        for (int i = crossover_point1; i <= crossover_point2; ++i){
            std::tie(child1.values[i], child2.values[i]) = crossover_rows(child1.values[i], child2.values[i]);
        }
    }

    return std::make_pair(child1, child2);
}

std::pair<std::vector<int>, std::vector<int>> CycleCrossover::crossover_rows(const std::vector<int> &row1, const std::vector<int> &row2){
    std::vector<int> child_row1(Nd, 0), child_row2(Nd, 0);
    std::vector<int> remaining(Nd);

    std::iota(remaining.begin(), remaining.end(), 1);
    int cycle = 0;

    while (std::count(child_row1.begin(), child_row1.end(), 0) > 0 && std::count(child_row2.begin(), child_row2.end(), 0) > 0){
        if (cycle % 2 == 0){
            int index = find_unused(row1, remaining);
            int start = row1[index];
            remaining.erase(std::remove(remaining.begin(), remaining.end(), row1[index]), remaining.end());
            child_row1[index] = row1[index];
            child_row2[index] = row2[index];
            int next = row2[index];

            while (next != start){
                index = find_value(row1, next);
                child_row1[index] = row1[index];
                remaining.erase(std::remove(remaining.begin(), remaining.end(), row1[index]), remaining.end());
                child_row2[index] = row2[index];
                next = row2[index];
            }
            cycle++;
        } else {
            int index = find_unused(row1, remaining);
            int start = row1[index];
            remaining.erase(std::remove(remaining.begin(), remaining.end(), row1[index]), remaining.end());
            child_row1[index] = row2[index];
            child_row2[index] = row1[index];
            int next = row2[index];

            while (next != start){
                index = find_value(row1, next);
                child_row1[index] = row2[index];
                remaining.erase(std::remove(remaining.begin(), remaining.end(), row1[index]), remaining.end());
                child_row2[index] = row1[index];
                next = row2[index];
            }
            cycle++;
        }
    }
    return std::make_pair(child_row1, child_row2);
}

int CycleCrossover::find_unused(const std::vector<int> &parent_row, const std::vector<int> &remaining){
    for (int i = 0; i < parent_row.size(); ++i){
        if (std::find(remaining.begin(), remaining.end(), parent_row[i]) != remaining.end()){
            return i;
        }
    }
    return -1;
}

int CycleCrossover::find_value(const std::vector<int> &parent_row, int value){
    auto it = std::find(parent_row.begin(), parent_row.end(), value);
    if (it != parent_row.end())    {
        return std::distance(parent_row.begin(), it);
    }
    return -1;
}

void Sudoku::load(const std::string &path){
    std::ifstream file(path);
    if (file.is_open()){
        std::vector<std::vector<int>> values(Nd, std::vector<int>(Nd));
        for (int i = 0; i < Nd; ++i){
            for (int j = 0; j < Nd; ++j){
                file >> values[i][j];
            }
        }
        given = Given(values);
    }
}

void Sudoku::save(const std::string &path, const Candidate &solution){
    std::ofstream file(path);
    if (file.is_open()){
        for (const auto &row : solution.values){
            for (const auto &val : row){
                file << val << " ";
            }
            file << "\n";
        }
    }
}

Candidate Sudoku::solve(){
    int Nc = 1000;
    int Ne = static_cast<int>(0.05 * Nc);
    int Ng = 1000;
    int Nm = 0;

    double phi = 0;
    double sigma = 1;
    // First case mutation
    double mutation_rate = 0.0;
    // Second case mutation
    //double mutation_rate = 0.06;
    // Third case mutation
    //double mutation_rate = 0.1;
    // Fourth case mutation
    //double mutation_rate = 0.5;
    // Fifth case mutation
    //double mutation_rate = 1.0;
    population.seed(Nc, given);
    
    int stale = 0;
    for (int generation = 0; generation < Ng; ++generation){
        std::cout << "Generation " << generation << std::endl;

        double best_fitness = 0.0;
        for (const auto &candidate : population.candidates){
            if (candidate.fitness == 1){
                std::cout << "Solution found at generation " << generation << "!" << std::endl;
                for (const auto &row : candidate.values){
                    for (const auto &val : row){
                        std::cout << val << " ";
                    }
                    std::cout << std::endl;
                }
                return candidate;
            }
            if (candidate.fitness > best_fitness){
                best_fitness = candidate.fitness;
            }
        }

        std::cout << "Best fitness: " << best_fitness << std::endl;

        Population next_population;
        population.sort();
        std::vector<Candidate> elites(population.candidates.begin(), population.candidates.begin() + Ne);

        for (int count = Ne; count < Nc; count += 2){
            Candidate parent1 = Tournament::compete(population.candidates);
            Candidate parent2 = Tournament::compete(population.candidates);

            auto children = CycleCrossover::crossover(parent1, parent2, 1.0);
            Candidate &child1 = children.first;
            Candidate &child2 = children.second;

            double old_fitness1 = child1.fitness;
            bool success1 = child1.mutate(mutation_rate, given);
            child1.update_fitness();
            if (success1){
                Nm++;
                if (child1.fitness > old_fitness1){
                    phi++;
                }
            }

            double old_fitness2 = child2.fitness;
            bool success2 = child2.mutate(mutation_rate, given);
            child2.update_fitness();
            if (success2){
                Nm++;
                if (child2.fitness > old_fitness2){
                    phi++;
                }
            }

            next_population.candidates.push_back(child1);
            next_population.candidates.push_back(child2);
        }

        for (const auto &elite : elites){
            next_population.candidates.push_back(elite);
        }

        population.candidates = next_population.candidates;
        population.update_fitness();

        if (Nm == 0){
            phi = 0;
        } else {
            phi /= Nm;
        }

        if (phi > 0.2){
            sigma /= 0.998;
        }
        else if (phi < 0.2){
            sigma *= 0.998;
        }

        static std::mt19937 rng(std::random_device{}());
        std::normal_distribution<> dist(0.0, sigma);
        mutation_rate = std::abs(dist(rng));
        Nm = 0;
        phi = 0;

        population.sort();
        if (population.candidates[0].fitness != population.candidates[1].fitness){
            stale = 0;
        } else {
            stale++;
        }

        if (stale >= 100){
            std::cout << "The population has gone stale. Re-seeding..." << std::endl;
            population.seed(Nc, given);
            stale = 0;
            sigma = 1;
            phi = 0;
            Nm = 0;
            // First case mutation
            mutation_rate = 0.0;
            // Second case mutation
            //double mutation_rate = 0.06;
            // Third case mutation
            //double mutation_rate = 0.1;
            // Fourth case mutation
            //double mutation_rate = 0.5;
            // Fifth case mutation
            //double mutation_rate = 1.0;
        }
    }

    std::cout << "No solution found." << std::endl;
    double best_fitness = 0.0;
    Candidate best_candidate;
    for (const auto &candidate : population.candidates){
            if (candidate.fitness > best_fitness){
                best_fitness = candidate.fitness;
                best_candidate = candidate;
            }
        }
    //return Candidate();
    return best_candidate;
}

int main(){
    Sudoku s;

    for (int i = 0; i < 10; i++){
        s.load("puzzle_mild.txt");
        Candidate solution = s.solve();
        if (solution.fitness == 1.0){
            s.save("solution.txt", solution);
        }else{
            s.save("partial_solution.txt",solution);
        }    
    }
    return 0;
}
