
# Substitution Cipher Decryption Tool (C++)

This project is a console application written in **C++** that automatically decrypts monoalphabetic substitution ciphers using a combination of:

- frequency analysis (monograms, bigrams, trigrams),
- dictionary matching,
- stochastic optimization algorithms (hill climbing, simulated annealing),
- and multi-threaded execution.

---

## 🔍 Features

- Language model with English frequency statistics and word dictionary
- Configurable decryption scoring based on monograms, bigrams, trigrams, and known words
- Supports:
  - Hill Climbing
  - Simulated Annealing
  - Optional Genetic Algorithm (template)
- Command-line options for input/output, verbosity, thread count
- Multithreading support for faster optimization

---

## 🚀 How to Use

### **Compile:**
```bash
g++ desifravimas_su_raktu.cpp -o decipher -std=c++17 -pthread
```

### **Run:**
```bash
./decipher -i encryptedtext.txt -o decryptedtext.txt -v -t 4
```

### **Options:**
- `-i`, `--input` – Encrypted input file (default: `encryptedtext.txt`)
- `-o`, `--output` – Output file for decrypted text (default: `decryptedtext.txt`)
- `-v`, `--verbose` – Enables debug and scoring output
- `-t`, `--threads` – Number of threads (defaults to system concurrency)
- `-h`, `--help` – Show usage help

---

## 📊 Algorithms Used

- **Frequency Analysis** (letter statistics compared to English)
- **Key optimization** via:
  - Hill Climbing
  - Simulated Annealing
  - (Extensible for Genetic Algorithms)

---

## 📁 Example Files

- `encryptedtext.txt` – Input file with ciphered text
- `decryptedtext.txt` – Output file with deciphered result

---

## 📌 Project Status

✅ Working prototype  
🧪 Supports multiple methods  
🔄 Potential for further improvement (GUI, multilingual models, training)

---

## 👤 Author

- **Your Name**
- GitHub: [https://github.com/yourusername](https://github.com/yourusername)

---

## 📜 License

MIT License – free to use and modify.
