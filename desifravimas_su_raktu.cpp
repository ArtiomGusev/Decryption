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
    struct KalbosModelis {
        std::string pavadinimas;
        std::unordered_map<std::string, double> monogramos;
        std::unordered_map<std::string, double> bigramos;
        std::unordered_map<std::string, double> trigramos;
        std::unordered_set<std::string> zodynas;
    };

    std::vector<KalbosModelis> kalbosModeliai;
    std::mt19937 rng;
    size_t maxTekstoDydis = 100000;
    size_t minDydis = 5000;
    int maxIteracijuBePagerinimo = 500;
    int maxVisoIteraciju = 50000;
    double pradineTemperatura = 10.0;
    double galutineTemperatura = 0.01;
    std::atomic<bool> stabdymas{false};

        //Konfigūracijos parametrai
    struct Konfiguracija {
        double monogramuSvoris = 1.0;
        double bigramuSvoris = 2.0;
        double trigramuSvoris = 4.0;
        double zodynoSvoris = 3.0;
        int gijuSk = std::thread::hardware_concurrency();
        bool naudotiSimuliuotaAtleidima = true;
        bool naudotiKopimuKala = true;
        bool naudotiGenetiniAlgoritma = false;
        int genetinesPopuliacijosDydis = 50;
        int tekstoTyrimoStrategija = 1; // 0=pradzia, 1=atsitiktine, 2=kelios
        bool adaptyvusSvoris = true;
        bool detalusIsvedimas = false;
    } konfig;

public:
    SubstitutionCipherSolver() : rng(std::chrono::steady_clock::now().time_since_epoch().count()) {
        inicializuotiKalbosModelius();
    }

    void inicializuotiKalbosModelius() {
        KalbosModelis anglu;
        anglu.pavadinimas = "Anglu";
        
        //Anglų kalbos raidžių dažnis
        anglu.monogramos = {
            {"E", 0.1202}, {"T", 0.0910}, {"A", 0.0812}, {"O", 0.0768}, {"I", 0.0731},
            {"N", 0.0695}, {"S", 0.0628}, {"R", 0.0602}, {"H", 0.0592}, {"D", 0.0432},
            {"L", 0.0398}, {"U", 0.0288}, {"C", 0.0271}, {"M", 0.0261}, {"F", 0.0230},
            {"Y", 0.0211}, {"W", 0.0209}, {"G", 0.0203}, {"P", 0.0182}, {"B", 0.0149},
            {"V", 0.0111}, {"K", 0.0069}, {"X", 0.0017}, {"Q", 0.0011}, {"J", 0.0010},
            {"Z", 0.0007}
        };

        anglu.bigramos = {
            {"TH", 2.71}, {"HE", 2.33}, {"IN", 2.03}, {"ER", 1.78}, {"AN", 1.61},
            {"RE", 1.41}, {"ND", 1.32}, {"AT", 1.21}, {"ON", 1.13}, {"NT", 1.12},
            {"HA", 1.08}, {"ES", 1.07}, {"ST", 1.05}, {"EN", 1.04}, {"ED", 1.02},
            {"TO", 0.99}, {"IT", 0.98}, {"OU", 0.94}, {"EA", 0.89}, {"HI", 0.87},
            {"IS", 0.86}, {"OR", 0.86}, {"TI", 0.85}, {"AS", 0.81}, {"TE", 0.79},
            {"ET", 0.76}, {"NG", 0.75}, {"OF", 0.74}, {"AL", 0.73}, {"DE", 0.73},
            {"SE", 0.73}, {"LE", 0.71}, {"SA", 0.69}, {"SI", 0.69}, {"AR", 0.68},
            {"VE", 0.68}, {"RA", 0.66}, {"LD", 0.65}, {"UR", 0.64}, {"NO", 0.63}
        };

        anglu.trigramos = {
            {"THE", 3.51}, {"AND", 1.59}, {"ING", 1.15}, {"HER", 0.82}, {"HAT", 0.65},
            {"HIS", 0.60}, {"THA", 0.59}, {"ERE", 0.55}, {"FOR", 0.55}, {"ENT", 0.53},
            {"ION", 0.51}, {"TER", 0.46}, {"WAS", 0.46}, {"YOU", 0.44}, {"ITH", 0.43},
            {"VER", 0.43}, {"ALL", 0.42}, {"WIT", 0.40}, {"THI", 0.39}, {"TIO", 0.38}
        };

        //Daugiausiai naudojami žodžiai
        anglu.zodynas = {
            "the", "be", "to", "of", "and", "a", "in", "that", "have", "i", "it", "for", "not", "on", "with",
            "he", "as", "you", "do", "at", "this", "but", "his", "by", "from", "they", "we", "say", "her", "she",
            "or", "an", "will", "my", "one", "all", "would", "there", "their", "what", "so", "up", "out", "if", "about",
            "who", "get", "which", "go", "me", "when", "make", "can", "like", "time", "no", "just", "him", "know", "take",
            "people", "into", "year", "your", "good", "some", "could", "them", "see", "other", "than", "then", "now", "look",
            "only", "come", "its", "over", "think", "also", "back", "after", "use", "two", "how", "our", "work", "first", "well",
            "way", "even", "new", "want", "because", "any", "these", "give", "day", "most", "us"
        };

        kalbosModeliai.push_back(anglu);
    }

    void nustatytiKonfiguracija(const Konfiguracija& naujaKonfig) {
        konfig = naujaKonfig;
    }

    Konfiguracija gautiKonfiguracija() const {
        return konfig;
    }

    void nustatytiMaxTekstoDydi(size_t dydis) {
        maxTekstoDydis = dydis;
    }

    void prasytiStabdymo() {
        stabdymas = true;
    }

    std::string desifruotiFaila(const std::string& ivestiesFailas, const std::string& isvestiesFailas, int kalbosModelioIndeksas = 0) {
        if (kalbosModelioIndeksas < 0 || kalbosModelioIndeksas >= kalbosModeliai.size()) {
            kalbosModelioIndeksas = 0;
        }

        try {
            std::string sifruotasTekstas = skaitytiFaila(ivestiesFailas);
            if (sifruotasTekstas.empty()) {
                std::cerr << "Klaida: Tuscias ivesties failas arba failas nerastas: " << ivestiesFailas << std::endl;
                return "";
            }

            auto pradzia = std::chrono::high_resolution_clock::now();

            if (konfig.detalusIsvedimas) {
                std::cout << "Teksto ilgis: " << sifruotasTekstas.size() << " simboliu" << std::endl;
                std::cout << "Naudojamas kalbos modelis: " << kalbosModeliai[kalbosModelioIndeksas].pavadinimas << std::endl;
                std::cout << "Pradedamas desifravimas..." << std::endl;
            }

            auto dazniuAnalize = atliktiDazniuAnalize(sifruotasTekstas);
            auto geriausiasRaktas = sukurtiPradiniRakta(dazniuAnalize, kalbosModelioIndeksas);
            auto geriausiasIvertinimas = ivertintiDesifravima(sifruotasTekstas, geriausiasRaktas, kalbosModelioIndeksas);
            
            if (konfig.detalusIsvedimas) {
                std::cout << "Pradinis rakto ivercinimas: " << geriausiasIvertinimas << std::endl;
            }

            geriausiasRaktas = optimizuotiRakta(sifruotasTekstas, geriausiasRaktas, kalbosModelioIndeksas);
            
             //Iššifravimas su geriausiu raktu
            std::string desifruotasTekstas = taikytiRakta(sifruotasTekstas, geriausiasRaktas);

            rasytiFaila(isvestiesFailas, desifruotasTekstas);

            auto pabaiga = std::chrono::high_resolution_clock::now();
            auto trukme = std::chrono::duration_cast<std::chrono::milliseconds>(pabaiga - pradzia).count();

            if (konfig.detalusIsvedimas) {
                std::cout << "Desifravimas baigtas per " << trukme / 1000.0 << " sekundziu" << std::endl;
                std::cout << "Desifruotas tekstas irasytas i: " << isvestiesFailas << std::endl;
                
                std::cout << "Galutinis rakto atvaizdavimas:" << std::endl;
                for (const auto& pora : geriausiasRaktas) {
                    std::cout << pora.first << " -> " << pora.second << "  ";
                }
                std::cout << std::endl;
            }

            return desifruotasTekstas;
        }
        catch (const std::exception& e) {
            std::cerr << "Klaida desifruojant: " << e.what() << std::endl;
            return "";
        }
    }

private:
    std::string skaitytiFaila(const std::string& failoVardas) {
        std::ifstream failas(failoVardas, std::ios::binary);
        if (!failas) {
            throw std::runtime_error("Nepavyko atidaryti failo: " + failoVardas);
        }

        //Gauti failo dydį
        failas.seekg(0, std::ios::end);
        size_t failoDydis = failas.tellg();
        failas.seekg(0, std::ios::beg);

        std::string turinys;
        turinys.resize(failoDydis);
        failas.read(&turinys[0], failoDydis);
        
        return turinys;
    }

    void rasytiFaila(const std::string& failoVardas, const std::string& turinys) {
        std::ofstream failas(failoVardas);
        if (!failas) {
            throw std::runtime_error("Nepavyko irasyti i faila: " + failoVardas);
        }
        failas << turinys;
    }

    std::unordered_map<char, double> atliktiDazniuAnalize(const std::string& tekstas) {
        std::unordered_map<char, double> dazniai;
        int visoSimboliu = 0;

        //Supaprastinti tekstą
        std::string tekstoTiriniui = tekstas;
        if (tekstas.size() > maxTekstoDydis) {
            tekstoTiriniui = tekstoTiriniui.substr(0, maxTekstoDydis);
        }

         //Nustatyti raidžių dažnumą
        for (char c : tekstoTiriniui) {
            if (std::isalpha(c)) {
                dazniai[std::toupper(c)]++;
                visoSimboliu++;
            }
        }

        if (visoSimboliu > 0) {
            for (auto& pora : dazniai) {
                pora.second /= visoSimboliu;
            }
        }

        return dazniai;
    }

    std::unordered_map<char, char> sukurtiPradiniRakta(
        const std::unordered_map<char, double>& dazniai, 
        int kalbosModelioIndeksas
    ) {
        std::vector<std::pair<char, double>> surikiuotiDazniai(dazniai.begin(), dazniai.end());
        std::sort(surikiuotiDazniai.begin(), surikiuotiDazniai.end(), 
            [](const auto& a, const auto& b) { return a.second > b.second; });

        std::vector<std::pair<std::string, double>> kalbosDazniai;
        for (const auto& pora : kalbosModeliai[kalbosModelioIndeksas].monogramos) {
            kalbosDazniai.push_back({pora.first, pora.second});
        }
        std::sort(kalbosDazniai.begin(), kalbosDazniai.end(), 
            [](const auto& a, const auto& b) { return a.second > b.second; });

        std::unordered_map<char, char> raktoZemelapis;
        size_t minDydis = std::min(surikiuotiDazniai.size(), kalbosDazniai.size());
        
        for (size_t i = 0; i < minDydis; ++i) {
            raktoZemelapis[surikiuotiDazniai[i].first] = kalbosDazniai[i].first[0];
        }

        for (char c = 'A'; c <= 'Z'; ++c) {
            if (raktoZemelapis.find(c) == raktoZemelapis.end()) {
                for (char p = 'A'; p <= 'Z'; ++p) {
                    bool naudotas = false;
                    for (const auto& pora : raktoZemelapis) {
                        if (pora.second == p) {
                            naudotas = true;
                            break;
                        }
                    }
                    if (!naudotas) {
                        raktoZemelapis[c] = p;
                        break;
                    }
                }
            }
        }

        return raktoZemelapis;
    }

    std::string taikytiRakta(const std::string& tekstas, const std::unordered_map<char, char>& raktoZemelapis) {
        std::string rezultatas;
        rezultatas.reserve(tekstas.size());

        for (char c : tekstas) {
            if (std::isalpha(c)) {
                char didzioji = std::toupper(c);
                auto it = raktoZemelapis.find(didzioji);
                if (it != raktoZemelapis.end()) {
                    rezultatas += std::islower(c) ? std::tolower(it->second) : it->second;
                } else {
                    rezultatas += c;
                }
            } else {
                rezultatas += c; 
            }
        }

        return rezultatas;
    }

     //Įvertinti dešifravimo kokybę
    double ivertintiDesifravima(
        const std::string& sifruotasTekstas,
        const std::unordered_map<char, char>& raktoZemelapis,
        int kalbosModelioIndeksas
    ) {
        std::string vertinimui = sifruotasTekstas;
        if (sifruotasTekstas.size() > maxTekstoDydis) {
            if (konfig.tekstoTyrimoStrategija == 0) {
                vertinimui = sifruotasTekstas.substr(0, maxTekstoDydis);
            } 
            else if (konfig.tekstoTyrimoStrategija == 1) {
                size_t pradziosPoz = rng() % (sifruotasTekstas.size() - maxTekstoDydis);
                vertinimui = sifruotasTekstas.substr(pradziosPoz, maxTekstoDydis);
            }
            else {
                size_t imtiesDydis = maxTekstoDydis / 3;
                std::string imtis1 = sifruotasTekstas.substr(0, imtiesDydis);
                size_t vidurys = sifruotasTekstas.size() / 2 - imtiesDydis / 2;
                std::string imtis2 = sifruotasTekstas.substr(vidurys, imtiesDydis);
                size_t pabaiga = sifruotasTekstas.size() - imtiesDydis;
                std::string imtis3 = sifruotasTekstas.substr(pabaiga, imtiesDydis);
                vertinimui = imtis1 + imtis2 + imtis3;
            }
        }

        std::string desifruotasTekstas = taikytiRakta(vertinimui, raktoZemelapis);

        double monogramuIvertis = 0.0;
        double bigramuIvertis = 0.0;
        double trigramuIvertis = 0.0;
        double zodynoIvertis = 0.0;

        std::unordered_map<std::string, int> monogramuSk;
        int visoMonogramu = 0;

        for (char c : desifruotasTekstas) {
            if (std::isalpha(c)) {
                std::string raide(1, std::toupper(c));
                monogramuSk[raide]++;
                visoMonogramu++;
            }
        }

        if (visoMonogramu > 0) {
            for (const auto& pora : monogramuSk) {
                double daznis = static_cast<double>(pora.second) / visoMonogramu;
                auto it = kalbosModeliai[kalbosModelioIndeksas].monogramos.find(pora.first);
                if (it != kalbosModeliai[kalbosModelioIndeksas].monogramos.end()) {
                    monogramuIvertis += (1.0 - std::abs(daznis - it->second));
                }
            }
            monogramuIvertis /= 26.0; 
        }

        std::unordered_map<std::string, int> bigramuSk;
        int visoBigramu = 0;

        for (size_t i = 1; i < desifruotasTekstas.size(); ++i) {
            if (std::isalpha(desifruotasTekstas[i-1]) && std::isalpha(desifruotasTekstas[i])) {
                std::string bigrama;
                bigrama += std::toupper(desifruotasTekstas[i-1]);
                bigrama += std::toupper(desifruotasTekstas[i]);
                bigramuSk[bigrama]++;
                visoBigramu++;
            }
        }

        if (visoBigramu > 0) {
            for (const auto& pora : bigramuSk) {
                auto it = kalbosModeliai[kalbosModelioIndeksas].bigramos.find(pora.first);
                if (it != kalbosModeliai[kalbosModelioIndeksas].bigramos.end()) {
                    bigramuIvertis += it->second * pora.second;
                }
            }
            bigramuIvertis /= visoBigramu;
        }

        std::unordered_map<std::string, int> trigramuSk;
        int visoTrigramu = 0;

        for (size_t i = 2; i < desifruotasTekstas.size(); ++i) {
            if (std::isalpha(desifruotasTekstas[i-2]) && std::isalpha(desifruotasTekstas[i-1]) && std::isalpha(desifruotasTekstas[i])) {
                std::string trigrama;
                trigrama += std::toupper(desifruotasTekstas[i-2]);
                trigrama += std::toupper(desifruotasTekstas[i-1]);
                trigrama += std::toupper(desifruotasTekstas[i]);
                trigramuSk[trigrama]++;
                visoTrigramu++;
            }
        }

        if (visoTrigramu > 0) {
            for (const auto& pora : trigramuSk) {
                auto it = kalbosModeliai[kalbosModelioIndeksas].trigramos.find(pora.first);
                if (it != kalbosModeliai[kalbosModelioIndeksas].trigramos.end()) {
                    trigramuIvertis += it->second * pora.second;
                }
            }
            trigramuIvertis /= visoTrigramu;
        }

        size_t poz = 0;
        int zodziuSk = 0;
        int rastuZodziu = 0;

        while (poz < desifruotasTekstas.size()) {
            while (poz < desifruotasTekstas.size() && !std::isalpha(desifruotasTekstas[poz])) {
                poz++;
            }
            
            if (poz >= desifruotasTekstas.size()) break;
            
            size_t zodzioPradzia = poz;
            while (poz < desifruotasTekstas.size() && std::isalpha(desifruotasTekstas[poz])) {
                poz++;
            }
            
            if (poz > zodzioPradzia) {
                zodziuSk++;
                std::string zodis = desifruotasTekstas.substr(zodzioPradzia, poz - zodzioPradzia);
                std::transform(zodis.begin(), zodis.end(), zodis.begin(), ::tolower);
                
                if (kalbosModeliai[kalbosModelioIndeksas].zodynas.find(zodis) != 
                    kalbosModeliai[kalbosModelioIndeksas].zodynas.end()) {
                    rastuZodziu++;
                }
            }
        }

        if (zodziuSk > 0) {
            zodynoIvertis = static_cast<double>(rastuZodziu) / zodziuSk;
        }

        double bendrasIvertis = 
            konfig.monogramuSvoris * monogramuIvertis + 
            konfig.bigramuSvoris * bigramuIvertis + 
            konfig.trigramuSvoris * trigramuIvertis + 
            konfig.zodynoSvoris * zodynoIvertis;

        bendrasIvertis = bendrasIvertis * 10.0 / (konfig.monogramuSvoris + konfig.bigramuSvoris + 
                                           konfig.trigramuSvoris + konfig.zodynoSvoris);

        return bendrasIvertis;
    }

    std::unordered_map<char, char> sukeistiDviRaites(
        const std::unordered_map<char, char>& raktoZemelapis
    ) {
        std::vector<char> raides;
        
        for (const auto& pora : raktoZemelapis) {
            raides.push_back(pora.second);
        }
        
        std::uniform_int_distribution<int> dist(0, raides.size() - 1);
        int poz1 = dist(rng);
        int poz2;
        do {
            poz2 = dist(rng);
        } while (poz1 == poz2);
        
        std::swap(raides[poz1], raides[poz2]);
        
        std::unordered_map<char, char> naujasRaktas;
        int i = 0;
        for (const auto& pora : raktoZemelapis) {
            naujasRaktas[pora.first] = raides[i++];
        }
        
        return naujasRaktas;
    }

    std::unordered_map<char, char> kopimuKala(
        const std::string& sifruotasTekstas,
        const std::unordered_map<char, char>& pradinisRaktas,
        int kalbosModelioIndeksas
    ) {
        auto geriausiasRaktas = pradinisRaktas;
        double geriausiasIvertis = ivertintiDesifravima(sifruotasTekstas, geriausiasRaktas, kalbosModelioIndeksas);
        
        int iteracijos = 0;
        int bePagerinimo = 0;
        
        while (iteracijos < maxVisoIteraciju && 
               bePagerinimo < maxIteracijuBePagerinimo &&
               !stabdymas) {
            
            auto naujasRaktas = sukeistiDviRaites(geriausiasRaktas);
            double naujasIvertis = ivertintiDesifravima(sifruotasTekstas, naujasRaktas, kalbosModelioIndeksas);
            
            if (naujasIvertis > geriausiasIvertis) {
                geriausiasRaktas = naujasRaktas;
                geriausiasIvertis = naujasIvertis;
                bePagerinimo = 0;
                
                if (konfig.detalusIsvedimas && iteracijos % 100 == 0) {
                    std::cout << "Pagerintas ivertinimas: " << geriausiasIvertis << " iteracijoje " << iteracijos << std::endl;
                }
            } else {
                bePagerinimo++;
            }
            
            iteracijos++;
        }
        
        if (konfig.detalusIsvedimas) {
            std::cout << "Kopimu kala baigta po " << iteracijos << " iteraciju" << std::endl;
            std::cout << "Galutinis ivertinimas: " << geriausiasIvertis << std::endl;
        }
        
        return geriausiasRaktas;
    }

    std::unordered_map<char, char> simuliatoriusAtleidimo(
        const std::string& sifruotasTekstas,
        const std::unordered_map<char, char>& pradinisRaktas,
        int kalbosModelioIndeksas
    ) {
        auto dabartinisRaktas = pradinisRaktas;
        double dabartinisIvertis = ivertintiDesifravima(sifruotasTekstas, dabartinisRaktas, kalbosModelioIndeksas);
        
        auto geriausiasRaktas = dabartinisRaktas;
        double geriausiasIvertis = dabartinisIvertis;
        
        int iteracijos = 0;
        double temperatura = pradineTemperatura;
        double atvesimoGreitis = std::pow(galutineTemperatura / pradineTemperatura, 1.0 / maxVisoIteraciju);
        
        while (iteracijos < maxVisoIteraciju && temperatura > galutineTemperatura && !stabdymas) {
            auto naujasRaktas = sukeistiDviRaites(dabartinisRaktas);
            double naujasIvertis = ivertintiDesifravima(sifruotasTekstas, naujasRaktas, kalbosModelioIndeksas);
            
            if (naujasIvertis > dabartinisIvertis) {
                dabartinisRaktas = naujasRaktas;
                dabartinisIvertis = naujasIvertis;
                
                if (dabartinisIvertis > geriausiasIvertis) {
                    geriausiasRaktas = dabartinisRaktas;
                    geriausiasIvertis = dabartinisIvertis;
                    
                    if (konfig.detalusIsvedimas && iteracijos % 100 == 0) {
                        std::cout << "Naujas geriausias ivertinimas: " << geriausiasIvertis 
                                  << " iteracijoje " << iteracijos 
                                  << ", temp: " << temperatura << std::endl;
                    }
                }
            } 
            else {
                double skirtumas = naujasIvertis - dabartinisIvertis;
                double priemimoTikimybe = std::exp(skirtumas / temperatura);
                std::uniform_real_distribution<double> dist(0.0, 1.0);
                
                if (dist(rng) < priemimoTikimybe) {
                    dabartinisRaktas = naujasRaktas;
                    dabartinisIvertis = naujasIvertis;
                }
            }
            
            temperatura *= atvesimoGreitis;
            iteracijos++;
        }
        
        if (konfig.detalusIsvedimas) {
            std::cout << "Simuliatorius atleidimo baigtas po " << iteracijos << " iteraciju" << std::endl;
            std::cout << "Galutinis geriausias ivertinimas: " << geriausiasIvertis << std::endl;
        }
        
        return geriausiasRaktas;
    }

    std::unordered_map<char, char> genetinisAlgoritmas(
        const std::string& sifruotasTekstas,
        const std::unordered_map<char, char>& pradinisRaktas,
        int kalbosModelioIndeksas
    ) {
        return pradinisRaktas;
    }

    std::unordered_map<char, char> optimizuotiRakta(
        const std::string& sifruotasTekstas,
        const std::unordered_map<char, char>& pradinisRaktas,
        int kalbosModelioIndeksas
    ) {
        std::vector<std::thread> gijos;
        std::vector<std::unordered_map<char, char>> rezultatai;
        std::vector<double> ivertinimai;
        std::mutex rezultatuMutex;
        
        rezultatai.resize(konfig.gijuSk);
        ivertinimai.resize(konfig.gijuSk, 0.0);
        
        for (int i = 0; i < konfig.gijuSk; ++i) {
            gijos.emplace_back([this, i, &sifruotasTekstas, &pradinisRaktas, 
                                 kalbosModelioIndeksas, &rezultatai, &ivertinimai, &rezultatuMutex]() {
                std::unordered_map<char, char> geriausiasRaktas;
                double geriausiasIvertis = 0.0;
                
                if (i % 3 == 0 && konfig.naudotiKopimuKala) {
                    geriausiasRaktas = kopimuKala(sifruotasTekstas, pradinisRaktas, kalbosModelioIndeksas);
                    geriausiasIvertis = ivertintiDesifravima(sifruotasTekstas, geriausiasRaktas, kalbosModelioIndeksas);
                } 
                else if (i % 3 == 1 && konfig.naudotiSimuliuotaAtleidima) {
                    geriausiasRaktas = simuliatoriusAtleidimo(sifruotasTekstas, pradinisRaktas, kalbosModelioIndeksas);
                    geriausiasIvertis = ivertintiDesifravima(sifruotasTekstas, geriausiasRaktas, kalbosModelioIndeksas);
                }
                else if (konfig.naudotiGenetiniAlgoritma) {
                    geriausiasRaktas = genetinisAlgoritmas(sifruotasTekstas, pradinisRaktas, kalbosModelioIndeksas);
                    geriausiasIvertis = ivertintiDesifravima(sifruotasTekstas, geriausiasRaktas, kalbosModelioIndeksas);
                }
                else {
                    geriausiasRaktas = kopimuKala(sifruotasTekstas, pradinisRaktas, kalbosModelioIndeksas);
                    geriausiasIvertis = ivertintiDesifravima(sifruotasTekstas, geriausiasRaktas, kalbosModelioIndeksas);
                }
                
                std::lock_guard<std::mutex> lock(rezultatuMutex);
                rezultatai[i] = geriausiasRaktas;
                ivertinimai[i] = geriausiasIvertis;
            });
        }
        
        for (auto& gija : gijos) {
            if (gija.joinable()) {
                gija.join();
            }
        }
        
        int geriausiasIndeksas = 0;
        double geriausiasIvertis = ivertinimai[0];
        
        for (int i = 1; i < ivertinimai.size(); ++i) {
            if (ivertinimai[i] > geriausiasIvertis) {
                geriausiasIvertis = ivertinimai[i];
                geriausiasIndeksas = i;
            }
        }
        
        if (konfig.detalusIsvedimas) {
            std::cout << "Geriausias optimizavimo rezultatas: " << geriausiasIvertis << std::endl;
        }
        
        return rezultatai[geriausiasIndeksas];
    }
};

int main(int argc, char* argv[]) {
    std::string ivestiesFailas = "uzsifruotastekstas.txt";
    std::string isvestiesFailas = "tikrastekstas.txt";
    bool detalusIsvedimas = false;
    int gijos = std::thread::hardware_concurrency();
    
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-i" || arg == "--input") {
            if (i + 1 < argc) {
                ivestiesFailas = argv[++i];
            }
        } else if (arg == "-o" || arg == "--output") {
            if (i + 1 < argc) {
                isvestiesFailas = argv[++i];
            }
        } else if (arg == "-v" || arg == "--verbose") {
            detalusIsvedimas = true;
        } else if (arg == "-t" || arg == "--threads") {
            if (i + 1 < argc) {
                gijos = std::stoi(argv[++i]);
            }
        } else if (arg == "-h" || arg == "--help") {
            std::cout << "Naudojimas: " << argv[0] << " [parinktys]" << std::endl;
            std::cout << "Parinktys:" << std::endl;
            std::cout << "  -i, --input FAILAS    Įvesties užšifruotas failas (numatytasis: uzsifruotastekstas.txt)" << std::endl;
            std::cout << "  -o, --output FAILAS   Išvesties dešifruotas failas (numatytasis: tikrastekstas.txt)" << std::endl;
            return 0;
        }
    }
    
    try {
        SubstitutionCipherSolver desifratorius;
        
        auto konfig = desifratorius.gautiKonfiguracija();
        konfig.detalusIsvedimas = detalusIsvedimas;
        konfig.gijuSk = gijos;
        desifratorius.nustatytiKonfiguracija(konfig);
        
        std::string desifruotasTekstas = desifratorius.desifruotiFaila(ivestiesFailas, isvestiesFailas);
        
        if (detalusIsvedimas) {
            size_t perziurosIlgis = std::min(desifruotasTekstas.size(), size_t(500));
            std::cout << "\nDesifruoto teksto perziura:\n" << desifruotasTekstas.substr(0, perziurosIlgis) << "..." << std::endl;
        } else {
            std::cout << "Desifravimas baigtas sekmingai. Rezultatas irasytas i " << isvestiesFailas << std::endl;
        }
        
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "Klaida: " << e.what() << std::endl;
        return 1;
    }
}
