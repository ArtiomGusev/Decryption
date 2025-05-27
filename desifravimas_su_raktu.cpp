#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <random>
#include <cctype>
#include <cmath>
#include <chrono>
#include <thread>
#include <mutex>
#include <atomic>
#include <functional>
#include <iomanip>
#include <sstream>
#include <memory>

class SubstitutionCipherSolver {
private:
    struct LanguageModel {
        std::string name;
        std::unordered_map<std::string, double> monograms;
        std::unordered_map<std::string, double> bigrams;
        std::unordered_map<std::string, double> trigrams;
        std::unordered_set<std::string> dictionary;
    };

    std::vector<LanguageModel> languageModels;
    std::mt19937 rng;
    size_t maxTextSize = 100000;
    size_t minSize = 5000;
    int maxIterationsWithoutImprovement = 500;
    int maxTotalIterations = 50000;
    double initialTemperature = 10.0;
    double finalTemperature = 0.01;
    std::atomic<bool> stopFlag{false};

    // Configuration parameters
    struct Configuration {
        double monogramWeight = 1.0;
        double bigramWeight = 2.0;
        double trigramWeight = 4.0;
        double dictionaryWeight = 3.0;
        int threadCount = std::thread::hardware_concurrency();
        bool useSimulatedAnnealing = true;
        bool useHillClimbing = true;
        bool useGeneticAlgorithm = false;
        int geneticPopulationSize = 50;
        int textSamplingStrategy = 1; // 0=beginning, 1=random, 2=multiple
        bool adaptiveWeight = true;
        bool verboseOutput = false;
    } config;

public:
    SubstitutionCipherSolver() : rng(std::chrono::steady_clock::now().time_since_epoch().count()) {
        initializeLanguageModels();
    }

    void initializeLanguageModels() {
        LanguageModel english;
        english.name = "English";
        
        // English letter frequencies
        english.monograms = {
            {"E", 0.1202}, {"T", 0.0910}, {"A", 0.0812}, {"O", 0.0768}, {"I", 0.0731},
            {"N", 0.0695}, {"S", 0.0628}, {"R", 0.0602}, {"H", 0.0592}, {"D", 0.0432},
            {"L", 0.0398}, {"U", 0.0288}, {"C", 0.0271}, {"M", 0.0261}, {"F", 0.0230},
            {"Y", 0.0211}, {"W", 0.0209}, {"G", 0.0203}, {"P", 0.0182}, {"B", 0.0149},
            {"V", 0.0111}, {"K", 0.0069}, {"X", 0.0017}, {"Q", 0.0011}, {"J", 0.0010},
            {"Z", 0.0007}
        };

        english.bigrams = {
            {"TH", 2.71}, {"HE", 2.33}, {"IN", 2.03}, {"ER", 1.78}, {"AN", 1.61},
            {"RE", 1.41}, {"ND", 1.32}, {"AT", 1.21}, {"ON", 1.13}, {"NT", 1.12},
            {"HA", 1.08}, {"ES", 1.07}, {"ST", 1.05}, {"EN", 1.04}, {"ED", 1.02},
            {"TO", 0.99}, {"IT", 0.98}, {"OU", 0.94}, {"EA", 0.89}, {"HI", 0.87},
            {"IS", 0.86}, {"OR", 0.86}, {"TI", 0.85}, {"AS", 0.81}, {"TE", 0.79},
            {"ET", 0.76}, {"NG", 0.75}, {"OF", 0.74}, {"AL", 0.73}, {"DE", 0.73},
            {"SE", 0.73}, {"LE", 0.71}, {"SA", 0.69}, {"SI", 0.69}, {"AR", 0.68},
            {"VE", 0.68}, {"RA", 0.66}, {"LD", 0.65}, {"UR", 0.64}, {"NO", 0.63}
        };

        english.trigrams = {
            {"THE", 3.51}, {"AND", 1.59}, {"ING", 1.15}, {"HER", 0.82}, {"HAT", 0.65},
            {"HIS", 0.60}, {"THA", 0.59}, {"ERE", 0.55}, {"FOR", 0.55}, {"ENT", 0.53},
            {"ION", 0.51}, {"TER", 0.46}, {"WAS", 0.46}, {"YOU", 0.44}, {"ITH", 0.43},
            {"VER", 0.43}, {"ALL", 0.42}, {"WIT", 0.40}, {"THI", 0.39}, {"TIO", 0.38}
        };

        // Most common words
        english.dictionary = {
            "the", "be", "to", "of", "and", "a", "in", "that", "have", "i", "it", "for", "not", "on", "with",
            "he", "as", "you", "do", "at", "this", "but", "his", "by", "from", "they", "we", "say", "her", "she",
            "or", "an", "will", "my", "one", "all", "would", "there", "their", "what", "so", "up", "out", "if", "about",
            "who", "get", "which", "go", "me", "when", "make", "can", "like", "time", "no", "just", "him", "know", "take",
            "people", "into", "year", "your", "good", "some", "could", "them", "see", "other", "than", "then", "now", "look",
            "only", "come", "its", "over", "think", "also", "back", "after", "use", "two", "how", "our", "work", "first", "well",
            "way", "even", "new", "want", "because", "any", "these", "give", "day", "most", "us"
        };

        languageModels.push_back(english);
    }

    void setConfiguration(const Configuration& newConfig) {
        config = newConfig;
    }

    Configuration getConfiguration() const {
        return config;
    }

    void setMaxTextSize(size_t size) {
        maxTextSize = size;
    }

    void requestStop() {
        stopFlag = true;
    }

    std::string decryptFile(const std::string& inputFile, const std::string& outputFile, int languageModelIndex = 0) {
        if (languageModelIndex < 0 || languageModelIndex >= languageModels.size()) {
            languageModelIndex = 0;
        }

        try {
            std::string encryptedText = readFile(inputFile);
            if (encryptedText.empty()) {
                std::cerr << "Error: Empty input file or file not found: " << inputFile << std::endl;
                return "";
            }

            auto startTime = std::chrono::high_resolution_clock::now();

            if (config.verboseOutput) {
                std::cout << "Text length: " << encryptedText.size() << " characters" << std::endl;
                std::cout << "Using language model: " << languageModels[languageModelIndex].name << std::endl;
                std::cout << "Starting decryption..." << std::endl;
            }

            auto frequencyAnalysis = performFrequencyAnalysis(encryptedText);
            auto bestKey = createInitialKey(frequencyAnalysis, languageModelIndex);
            auto bestScore = evaluateDecryption(encryptedText, bestKey, languageModelIndex);
            
            if (config.verboseOutput) {
                std::cout << "Initial key score: " << bestScore << std::endl;
            }

            bestKey = optimizeKey(encryptedText, bestKey, languageModelIndex);
            
            // Decrypt with the best key
            std::string decryptedText = applyKey(encryptedText, bestKey);

            writeFile(outputFile, decryptedText);

            auto endTime = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();

            if (config.verboseOutput) {
                std::cout << "Decryption completed in " << duration / 1000.0 << " seconds" << std::endl;
                std::cout << "Decrypted text written to: " << outputFile << std::endl;
                
                std::cout << "Final key mapping:" << std::endl;
                for (const auto& pair : bestKey) {
                    std::cout << pair.first << " -> " << pair.second << "  ";
                }
                std::cout << std::endl;
            }

            return decryptedText;
        }
        catch (const std::exception& e) {
            std::cerr << "Error decrypting: " << e.what() << std::endl;
            return "";
        }
    }

private:
    std::string readFile(const std::string& filename) {
        std::ifstream file(filename, std::ios::binary);
        if (!file) {
            throw std::runtime_error("Failed to open file: " + filename);
        }

        // Get file size
        file.seekg(0, std::ios::end);
        size_t fileSize = file.tellg();
        file.seekg(0, std::ios::beg);

        std::string content;
        content.resize(fileSize);
        file.read(&content[0], fileSize);
        
        return content;
    }

    void writeFile(const std::string& filename, const std::string& content) {
        std::ofstream file(filename);
        if (!file) {
            throw std::runtime_error("Failed to write to file: " + filename);
        }
        file << content;
    }

    std::unordered_map<char, double> performFrequencyAnalysis(const std::string& text) {
        std::unordered_map<char, double> frequencies;
        int totalChars = 0;

        // Simplify text
        std::string textToAnalyze = text;
        if (text.size() > maxTextSize) {
            textToAnalyze = textToAnalyze.substr(0, maxTextSize);
        }

        // Calculate letter frequencies
        for (char c : textToAnalyze) {
            if (std::isalpha(c)) {
                frequencies[std::toupper(c)]++;
                totalChars++;
            }
        }

        if (totalChars > 0) {
            for (auto& pair : frequencies) {
                pair.second /= totalChars;
            }
        }

        return frequencies;
    }

    std::unordered_map<char, char> createInitialKey(
        const std::unordered_map<char, double>& frequencies, 
        int languageModelIndex
    ) {
        std::vector<std::pair<char, double>> sortedFrequencies(frequencies.begin(), frequencies.end());
        std::sort(sortedFrequencies.begin(), sortedFrequencies.end(), 
            [](const auto& a, const auto& b) { return a.second > b.second; });

        std::vector<std::pair<std::string, double>> languageFrequencies;
        for (const auto& pair : languageModels[languageModelIndex].monograms) {
            languageFrequencies.push_back({pair.first, pair.second});
        }
        std::sort(languageFrequencies.begin(), languageFrequencies.end(), 
            [](const auto& a, const auto& b) { return a.second > b.second; });

        std::unordered_map<char, char> keyMapping;
        size_t minSize = std::min(sortedFrequencies.size(), languageFrequencies.size());
        
        for (size_t i = 0; i < minSize; ++i) {
            keyMapping[sortedFrequencies[i].first] = languageFrequencies[i].first[0];
        }

        // Fill remaining letters with unused mappings
        for (char c = 'A'; c <= 'Z'; ++c) {
            if (keyMapping.find(c) == keyMapping.end()) {
                for (char p = 'A'; p <= 'Z'; ++p) {
                    bool used = false;
                    for (const auto& pair : keyMapping) {
                        if (pair.second == p) {
                            used = true;
                            break;
                        }
                    }
                    if (!used) {
                        keyMapping[c] = p;
                        break;
                    }
                }
            }
        }

        return keyMapping;
    }

    std::string applyKey(const std::string& text, const std::unordered_map<char, char>& keyMapping) {
        std::string result;
        result.reserve(text.size());

        for (char c : text) {
            if (std::isalpha(c)) {
                char upper = std::toupper(c);
                auto it = keyMapping.find(upper);
                if (it != keyMapping.end()) {
                    result += std::islower(c) ? std::tolower(it->second) : it->second;
                } else {
                    result += c;
                }
            } else {
                result += c; 
            }
        }

        return result;
    }

    // Evaluate decryption quality
    double evaluateDecryption(
        const std::string& encryptedText,
        const std::unordered_map<char, char>& keyMapping,
        int languageModelIndex
    ) {
        std::string textToEvaluate = encryptedText;
        if (encryptedText.size() > maxTextSize) {
            if (config.textSamplingStrategy == 0) {
                textToEvaluate = encryptedText.substr(0, maxTextSize);
            } 
            else if (config.textSamplingStrategy == 1) {
                size_t startPos = rng() % (encryptedText.size() - maxTextSize);
                textToEvaluate = encryptedText.substr(startPos, maxTextSize);
            }
            else {
                size_t sampleSize = maxTextSize / 3;
                std::string sample1 = encryptedText.substr(0, sampleSize);
                size_t middle = encryptedText.size() / 2 - sampleSize / 2;
                std::string sample2 = encryptedText.substr(middle, sampleSize);
                size_t end = encryptedText.size() - sampleSize;
                std::string sample3 = encryptedText.substr(end, sampleSize);
                textToEvaluate = sample1 + sample2 + sample3;
            }
        }

        std::string decryptedText = applyKey(textToEvaluate, keyMapping);

        double monogramScore = 0.0;
        double bigramScore = 0.0;
        double trigramScore = 0.0;
        double dictionaryScore = 0.0;

        std::unordered_map<std::string, int> monogramCounts;
        int totalMonograms = 0;

        for (char c : decryptedText) {
            if (std::isalpha(c)) {
                std::string letter(1, std::toupper(c));
                monogramCounts[letter]++;
                totalMonograms++;
            }
        }

        if (totalMonograms > 0) {
            for (const auto& pair : monogramCounts) {
                double frequency = static_cast<double>(pair.second) / totalMonograms;
                auto it = languageModels[languageModelIndex].monograms.find(pair.first);
                if (it != languageModels[languageModelIndex].monograms.end()) {
                    monogramScore += (1.0 - std::abs(frequency - it->second));
                }
            }
            monogramScore /= 26.0; 
        }

        std::unordered_map<std::string, int> bigramCounts;
        int totalBigrams = 0;

        for (size_t i = 1; i < decryptedText.size(); ++i) {
            if (std::isalpha(decryptedText[i-1]) && std::isalpha(decryptedText[i])) {
                std::string bigram;
                bigram += std::toupper(decryptedText[i-1]);
                bigram += std::toupper(decryptedText[i]);
                bigramCounts[bigram]++;
                totalBigrams++;
            }
        }

        if (totalBigrams > 0) {
            for (const auto& pair : bigramCounts) {
                auto it = languageModels[languageModelIndex].bigrams.find(pair.first);
                if (it != languageModels[languageModelIndex].bigrams.end()) {
                    bigramScore += it->second * pair.second;
                }
            }
            bigramScore /= totalBigrams;
        }

        std::unordered_map<std::string, int> trigramCounts;
        int totalTrigrams = 0;

        for (size_t i = 2; i < decryptedText.size(); ++i) {
            if (std::isalpha(decryptedText[i-2]) && std::isalpha(decryptedText[i-1]) && std::isalpha(decryptedText[i])) {
                std::string trigram;
                trigram += std::toupper(decryptedText[i-2]);
                trigram += std::toupper(decryptedText[i-1]);
                trigram += std::toupper(decryptedText[i]);
                trigramCounts[trigram]++;
                totalTrigrams++;
            }
        }

        if (totalTrigrams > 0) {
            for (const auto& pair : trigramCounts) {
                auto it = languageModels[languageModelIndex].trigrams.find(pair.first);
                if (it != languageModels[languageModelIndex].trigrams.end()) {
                    trigramScore += it->second * pair.second;
                }
            }
            trigramScore /= totalTrigrams;
        }

        size_t pos = 0;
        int wordCount = 0;
        int matchedWords = 0;

        while (pos < decryptedText.size()) {
            while (pos < decryptedText.size() && !std::isalpha(decryptedText[pos])) {
                pos++;
            }
            
            if (pos >= decryptedText.size()) break;
            
            size_t wordStart = pos;
            while (pos < decryptedText.size() && std::isalpha(decryptedText[pos])) {
                pos++;
            }
            
            if (pos > wordStart) {
                wordCount++;
                std::string word = decryptedText.substr(wordStart, pos - wordStart);
                std::transform(word.begin(), word.end(), word.begin(), ::tolower);
                
                if (languageModels[languageModelIndex].dictionary.find(word) != 
                    languageModels[languageModelIndex].dictionary.end()) {
                    matchedWords++;
                }
            }
        }

        if (wordCount > 0) {
            dictionaryScore = static_cast<double>(matchedWords) / wordCount;
        }

        double totalScore = 
            config.monogramWeight * monogramScore + 
            config.bigramWeight * bigramScore + 
            config.trigramWeight * trigramScore + 
            config.dictionaryWeight * dictionaryScore;

        totalScore = totalScore * 10.0 / (config.monogramWeight + config.bigramWeight + 
                                           config.trigramWeight + config.dictionaryWeight);

        return totalScore;
    }

    std::unordered_map<char, char> swapTwoLetters(
        const std::unordered_map<char, char>& keyMapping
    ) {
        std::vector<char> letters;
        
        for (const auto& pair : keyMapping) {
            letters.push_back(pair.second);
        }
        
        std::uniform_int_distribution<int> dist(0, letters.size() - 1);
        int pos1 = dist(rng);
        int pos2;
        do {
            pos2 = dist(rng);
        } while (pos1 == pos2);
        
        std::swap(letters[pos1], letters[pos2]);
        
        std::unordered_map<char, char> newKey;
        int i = 0;
        for (const auto& pair : keyMapping) {
            newKey[pair.first] = letters[i++];
        }
        
        return newKey;
    }

    std::unordered_map<char, char> hillClimbing(
        const std::string& encryptedText,
        const std::unordered_map<char, char>& initialKey,
        int languageModelIndex
    ) {
        auto bestKey = initialKey;
        double bestScore = evaluateDecryption(encryptedText, bestKey, languageModelIndex);
        
        int iterations = 0;
        int noImprovementCount = 0;
        
        while (iterations < maxTotalIterations && 
               noImprovementCount < maxIterationsWithoutImprovement &&
               !stopFlag) {
            
            auto newKey = swapTwoLetters(bestKey);
            double newScore = evaluateDecryption(encryptedText, newKey, languageModelIndex);
            
            if (newScore > bestScore) {
                bestKey = newKey;
                bestScore = newScore;
                noImprovementCount = 0;
                
                if (config.verboseOutput && iterations % 100 == 0) {
                    std::cout << "Improved score: " << bestScore << " at iteration " << iterations << std::endl;
                }
            } else {
                noImprovementCount++;
            }
            
            iterations++;
        }
        
        if (config.verboseOutput) {
            std::cout << "Hill climbing completed after " << iterations << " iterations" << std::endl;
            std::cout << "Final score: " << bestScore << std::endl;
        }
        
        return bestKey;
    }

    std::unordered_map<char, char> simulatedAnnealing(
        const std::string& encryptedText,
        const std::unordered_map<char, char>& initialKey,
        int languageModelIndex
    ) {
        auto currentKey = initialKey;
        double currentScore = evaluateDecryption(encryptedText, currentKey, languageModelIndex);
        
        auto bestKey = currentKey;
        double bestScore = currentScore;
        
        int iterations = 0;
        double temperature = initialTemperature;
        double coolingRate = std::pow(finalTemperature / initialTemperature, 1.0 / maxTotalIterations);
        
        while (iterations < maxTotalIterations && temperature > finalTemperature && !stopFlag) {
            auto newKey = swapTwoLetters(currentKey);
            double newScore = evaluateDecryption(encryptedText, newKey, languageModelIndex);
            
            if (newScore > currentScore) {
                currentKey = newKey;
                currentScore = newScore;
                
                if (currentScore > bestScore) {
                    bestKey = currentKey;
                    bestScore = currentScore;
                    
                    if (config.verboseOutput && iterations % 100 == 0) {
                        std::cout << "New best score: " << bestScore 
                                  << " at iteration " << iterations 
                                  << ", temp: " << temperature << std::endl;
                    }
                }
            } 
            else {
                double difference = newScore - currentScore;
                double acceptanceProbability = std::exp(difference / temperature);
                std::uniform_real_distribution<double> dist(0.0, 1.0);
                
                if (dist(rng) < acceptanceProbability) {
                    currentKey = newKey;
                    currentScore = newScore;
                }
            }
            
            temperature *= coolingRate;
            iterations++;
        }
        
        if (config.verboseOutput) {
            std::cout << "Simulated annealing completed after " << iterations << " iterations" << std::endl;
            std::cout << "Final best score: " << bestScore << std::endl;
        }
        
        return bestKey;
    }

    std::unordered_map<char, char> geneticAlgorithm(
        const std::string& encryptedText,
        const std::unordered_map<char, char>& initialKey,
        int languageModelIndex
    ) {
        return initialKey;
    }

    std::unordered_map<char, char> optimizeKey(
        const std::string& encryptedText,
        const std::unordered_map<char, char>& initialKey,
        int languageModelIndex
    ) {
        std::vector<std::thread> threads;
        std::vector<std::unordered_map<char, char>> results;
        std::vector<double> scores;
        std::mutex resultsMutex;
        
        results.resize(config.threadCount);
        scores.resize(config.threadCount, 0.0);
        
        for (int i = 0; i < config.threadCount; ++i) {
            threads.emplace_back([this, i, &encryptedText, &initialKey, 
                                 languageModelIndex, &results, &scores, &resultsMutex]() {
                std::unordered_map<char, char> bestKey;
                double bestScore = 0.0;
                
                if (i % 3 == 0 && config.useHillClimbing) {
                    bestKey = hillClimbing(encryptedText, initialKey, languageModelIndex);
                    bestScore = evaluateDecryption(encryptedText, bestKey, languageModelIndex);
                } 
                else if (i % 3 == 1 && config.useSimulatedAnnealing) {
                    bestKey = simulatedAnnealing(encryptedText, initialKey, languageModelIndex);
                    bestScore = evaluateDecryption(encryptedText, bestKey, languageModelIndex);
                }
                else if (config.useGeneticAlgorithm) {
                    bestKey = geneticAlgorithm(encryptedText, initialKey, languageModelIndex);
                    bestScore = evaluateDecryption(encryptedText, bestKey, languageModelIndex);
                }
                else {
                    bestKey = hillClimbing(encryptedText, initialKey, languageModelIndex);
                    bestScore = evaluateDecryption(encryptedText, bestKey, languageModelIndex);
                }
                
                std::lock_guard<std::mutex> lock(resultsMutex);
                results[i] = bestKey;
                scores[i] = bestScore;
            });
        }
        
        for (auto& thread : threads) {
            if (thread.joinable()) {
                thread.join();
            }
        }
        
        int bestIndex = 0;
        double bestScore = scores[0];
        
        for (int i = 1; i < scores.size(); ++i) {
            if (scores[i] > bestScore) {
                bestScore = scores[i];
                bestIndex = i;
            }
        }
        
        if (config.verboseOutput) {
            std::cout << "Best optimization result: " << bestScore << std::endl;
        }
        
        return results[bestIndex];
    }
};

int main(int argc, char* argv[]) {
    std::string inputFile = "encryptedtext.txt";
    std::string outputFile = "decryptedtext.txt";
    bool verboseOutput = false;
    int threads = std::thread::hardware_concurrency();
    
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-i" || arg == "--input") {
            if (i + 1 < argc) {
                inputFile = argv[++i];
            }
        } else if (arg == "-o" || arg == "--output") {
            if (i + 1 < argc) {
                outputFile = argv[++i];
            }
        } else if (arg == "-v" || arg == "--verbose") {
            verboseOutput = true;
        } else if (arg == "-t" || arg == "--threads") {
            if (i + 1 < argc) {
                threads = std::stoi(argv[++i]);
            }
        } else if (arg == "-h" || arg == "--help") {
            std::cout << "Usage: " << argv[0] << " [options]" << std::endl;
            std::cout << "Options:" << std::endl;
            std::cout << "  -i, --input FILE    Input encrypted file (default: encryptedtext.txt)" << std::endl;
            std::cout << "  -o, --output FILE   Output decrypted file (default: decryptedtext.txt)" << std::endl;
            std::cout << "  -v, --verbose       Enable verbose output" << std::endl;
            std::cout << "  -t, --threads N     Set number of threads to use" << std::endl;
            std::cout << "  -h, --help          Show this help message" << std::endl;
            return 0;
        }
    }
    
    try {
        SubstitutionCipherSolver decryptor;
        
        auto config = decryptor.getConfiguration();
        config.verboseOutput = verboseOutput;
        config.threadCount = threads;
        decryptor.setConfiguration(config);
        
        std::string decryptedText = decryptor.decryptFile(inputFile, outputFile);
        
        if (verboseOutput) {
            size_t previewLength = std::min(decryptedText.size(), size_t(500));
            std::cout << "\nDecrypted text preview:\n" << decryptedText.substr(0, previewLength) << "..." << std::endl;
        } else {
            std::cout << "Decryption completed successfully. Result written to " << outputFile << std::endl;
        }
        
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}