# Szczegółowe wyjaśnienie działania programu optymalizacji (Konspekt K3)

## Spis treści
1. [Wprowadzenie](#wprowadzenie)
2. [Struktura projektu](#struktura-projektu)
3. [Poprzednie konspekty (K0, K1, K2)](#poprzednie-konspekty)
4. [Konspekt K3 - Optymalizacja z ograniczeniami](#konspekt-k3)
5. [Algorytm sympleksu Neldera-Meada](#algorytm-nelder-mead)
6. [Algorytmy funkcji kary](#algorytmy-funkcji-kary)
7. [Problem testowy vs rzeczywisty](#problem-testowy-vs-rzeczywisty)
8. [Analiza istniejącej implementacji](#analiza-implementacji)

---

## Wprowadzenie

Program stanowi implementację różnych metod optymalizacji bezgradientowych do wyznaczania minimum funkcji celu. W konspekcie K3 głównym celem jest **optymalizacja z ograniczeniami** wykorzystująca metodę sympleks Neldera-Meada oraz **funkcje kary** (zewnętrzną i wewnętrzną).

---

## Struktura projektu

```
OPT/
├── main.cpp          # Główny plik z funkcjami lab0-lab6
├── opt_alg.cpp       # Implementacje algorytmów optymalizacji
├── opt_alg.h         # Deklaracje algorytmów
├── solution.cpp      # Klasa reprezentująca rozwiązanie
├── solution.h        # Deklaracja klasy solution
├── user_funs.cpp     # Funkcje celu (testowe i rzeczywiste)
├── user_funs.h       # Deklaracje funkcji celu
├── matrix.cpp        # Operacje na macierzach
├── matrix.h          # Deklaracje klasy matrix
├── ode_solver.cpp    # Solver równań różniczkowych (RK4)
├── ode_solver.h      # Deklaracja solvera
└── K3.pdf            # Konspekt do ćwiczenia 3
```

---

## Poprzednie konspekty

### Lab 0 - Monte Carlo i wahadło
```cpp
void lab0()
{
    // 1. Testowa funkcja celu - szukanie minimum f(x) = (x₁-a₁)² + (x₂-a₂)²
    solution opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);
    
    // 2. Problem rzeczywisty - optymalizacja momentu siły wahadła
    // Celem jest znalezienie takiego momentu M, który spowoduje
    // maksymalne wychylenie wahadła równe zadanej wartości teta_opt
    opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);
}
```

**Algorytm Monte Carlo (MC)**:
- Losuje punkty w przestrzeni rozwiązań [lb, ub]
- Oblicza wartość funkcji celu dla każdego punktu
- Zatrzymuje się gdy y < epsilon lub przekroczy Nmax wywołań

### Lab 1 - Metody jednowymiarowe
```cpp
void lab1(int aN)
{
    // Ekspansja - wyznaczanie przedziału zawierającego minimum
    p = expansion(ff1R, x0, d, alpha[aN], Nmax, cel);
    
    // Fibonacci - podział złotego ciągu
    opt = fib(ff1R, p[0], p[1], epsilon, cel);
    
    // Lagrange - interpolacja kwadratowa
    opt = lag(ff1R, ld, ud, epsilon, gamma, Nmax, cel, NAN);
}
```

**Problem rzeczywisty K1**: System dwóch zbiorników cieczy. Celem jest znalezienie takiej średnicy otworu D_A, aby maksymalna temperatura w zbiorniku B wyniosła 50°C.

### Lab 2 - Metody wielowymiarowe bez gradientu
```cpp
void lab2()
{
    // Metoda Hooke'a-Jeevesa (pattern search)
    optHJ = HJ(ff2R, k0, s, alphaHJ, epsilon, Nmax);
    
    // Metoda Rosenbrocka
    optRos = Rosen(ff2R, k0, s0, alphaRos, beta, epsilon, Nmax);
}
```

**Problem rzeczywisty K2**: Regulator PD dla ramienia robota. Celem jest znalezienie optymalnych współczynników k₁, k₂ regulatora.

---

## Konspekt K3 - Optymalizacja z ograniczeniami

### Testowa funkcja celu

Funkcja:
```
f(x₁, x₂) = sin(π√((x₁/π)² + (x₂/π)²)) / (π√((x₁/π)² + (x₂/π)²))
```

Implementacja w kodzie (`user_funs.cpp`):
```cpp
matrix ff3T(matrix x, matrix ud1, matrix ud2) {
    double x1 = m2d(x(0));
    double x2 = m2d(x(1));
    double common = M_PI * sqrt(pow(x1/M_PI, 2) + pow(x2/M_PI, 2));
    double t1 = sin(common);
    return matrix(t1 / common);
}
```

### Ograniczenia

Dla testowej funkcji celu:
1. **g₁(x₁) = -x₁ + 1 ≤ 0** → x₁ ≥ 1
2. **g₂(x₂) = -x₂ + 1 ≤ 0** → x₂ ≥ 1  
3. **g₃(x₁, x₂) = √(x₁² + x₂²) - a ≤ 0** → odległość od początku ≤ a

Gdzie `a` przyjmuje wartości: 4, 4.4934, lub 5.

### Problem rzeczywisty K3

**Scenariusz**: Piłka o masie m=600g i promieniu r=12cm spada z wysokości y₀=100m. Ma początkową prędkość poziomą v₀ₓ oraz rotację ω. Występuje efekt Magnusa.

**Cel optymalizacji**: Znaleźć v₀ₓ ∈ [-10, 10] m/s i ω ∈ [-10, 10] rad/s, które maksymalizują xₑₙd (odległość poziomą w momencie upadku).

**Ograniczenie**: Piłka musi minąć punkt (5, 50) w odległości nie większej niż 2m (dla y=50: x ∈ [3, 7]).

---

## Algorytm sympleksu Neldera-Meada

### Idea algorytmu

Metoda Nelder-Mead operuje na **sympleksie** - figurze geometrycznej w n-wymiarowej przestrzeni. Dla 2D jest to trójkąt, dla 3D czworościan.

### Pseudokod
```
1. Utwórz sympleks początkowy (n+1 wierzchołków)
   p₀ = x⁽⁰⁾ (punkt startowy)
   pᵢ = p₀ + s·eᵢ dla i = 1,...,n

2. Powtarzaj:
   a) Oblicz f(p₀), f(p₁), ..., f(pₙ)
   b) Znajdź p_min (najlepszy) i p_max (najgorszy)
   c) Oblicz centroid (środek ciężkości bez p_max):
      p̄ = (Σᵢ≠max pᵢ) / n
   
   d) ODBICIE: p_odb = p̄ + α(p̄ - p_max)
   
   e) Jeśli f(p_odb) < f(p_min):
      - EKSPANSJA: pₑ = p̄ + γ(p_odb - p̄)
      - Jeśli f(pₑ) < f(p_odb): p_max = pₑ
      - W przeciwnym razie: p_max = p_odb
      
   f) Jeśli f(p_min) ≤ f(p_odb) < f(p_max):
      - p_max = p_odb (akceptuj odbicie)
      
   g) Jeśli f(p_odb) ≥ f(p_max):
      - ZAWĘŻENIE: pᵤ = p̄ + β(p_max - p̄)
      - Jeśli f(pᵤ) < f(p_max): p_max = pᵤ
      - W przeciwnym razie: REDUKCJA wszystkich punktów w kierunku p_min
        pᵢ = δ(pᵢ + p_min) dla i ≠ min

3. Aż max||p_min - pᵢ|| < ε
```

### Parametry standardowe
- **α = 1** - współczynnik odbicia
- **β = 0.5** - współczynnik zawężenia
- **γ = 2** - współczynnik ekspansji
- **δ = 0.5** - współczynnik redukcji

### Implementacja w kodzie (`opt_alg.cpp`):

```cpp
solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, 
                double alpha, double beta, double gamma, double delta, 
                double epsilon, int Nmax, matrix ud1, matrix ud2)
{
    solution Xopt;
    matrix p[3] = {matrix(2,1), matrix(2,1), matrix(2,1)};
    p[0] = x0;  // Pierwszy wierzchołek = punkt startowy
    
    // Tworzenie sympleksu początkowego (n+1 wierzchołków dla n wymiarów)
    // p[0] = x0 (już ustawione)
    // p[1] = x0 + s*e1 (e1 = [1,0]ᵀ)
    // p[2] = x0 + s*e2 (e2 = [0,1]ᵀ)
    for (int i = 0; i < 2; i++) {
        p[i+1] = x0 + s * ident_mat(2)[i];  // UWAGA: p[i+1], nie p[i]!
    }
    
    // ... (główna pętla algorytmu)
}
```

---

## Algorytmy funkcji kary

### Idea metody funkcji kary

Metoda funkcji kary przekształca problem optymalizacji z ograniczeniami w problem bez ograniczeń poprzez dodanie do funkcji celu **kary** za naruszenie ograniczeń.

**Ogólna postać**:
```
F(x) = f(x) + c · S(x)
```

gdzie:
- `f(x)` - oryginalna funkcja celu
- `c` - współczynnik kary (zmieniany iteracyjnie)
- `S(x)` - funkcja kary

### Zewnętrzna funkcja kary (Exterior Penalty Method)

**Formuła**:
```
S(x₁, x₂) = Σᵢ (max(0, gᵢ(x₁, x₂)))²
```

**Jak działa**:
1. Jeśli ograniczenie jest spełnione (gᵢ ≤ 0): kara = 0
2. Jeśli ograniczenie jest naruszone (gᵢ > 0): kara = gᵢ²

**Przykład dla K3**:
```cpp
// Dla ograniczeń:
// g₁ = -x₁ + 1 ≤ 0  (x₁ ≥ 1)
// g₂ = -x₂ + 1 ≤ 0  (x₂ ≥ 1)
// g₃ = √(x₁² + x₂²) - a ≤ 0

double S_zewn(double x1, double x2, double a) {
    double g1 = -x1 + 1;  // ≤ 0 gdy x1 ≥ 1
    double g2 = -x2 + 1;  // ≤ 0 gdy x2 ≥ 1
    double g3 = sqrt(x1*x1 + x2*x2) - a;  // ≤ 0 gdy wewnątrz koła
    
    return pow(max(0.0, g1), 2) + 
           pow(max(0.0, g2), 2) + 
           pow(max(0.0, g3), 2);
}
```

**Charakterystyka**:
- Punkt startowy może być **poza** obszarem dopuszczalnym
- c zaczyna od małej wartości i **rośnie** (np. c(i+1) = α·c(i), α > 1)
- Rozwiązanie zbliża się do brzegu od zewnątrz

### Wewnętrzna funkcja kary (Interior Penalty / Barrier Method)

**Formuła**:
```
S(x₁, x₂) = -Σᵢ (1 / gᵢ(x₁, x₂))
```

**Jak działa**:
1. Gdy gᵢ → 0⁻ (zbliżamy się do granicy): kara → +∞
2. Gdy gᵢ << 0 (daleko od granicy): kara → 0

**Przykład dla K3**:
```cpp
double S_wewn(double x1, double x2, double a) {
    double g1 = -x1 + 1;
    double g2 = -x2 + 1;
    double g3 = sqrt(x1*x1 + x2*x2) - a;
    
    // Wszystkie g muszą być < 0 (wewnątrz obszaru)
    if (g1 >= 0 || g2 >= 0 || g3 >= 0) {
        return INFINITY;  // Punkt niedopuszczalny
    }
    
    return -1.0/g1 - 1.0/g2 - 1.0/g3;
}
```

**Charakterystyka**:
- Punkt startowy **musi być** wewnątrz obszaru dopuszczalnego
- c zaczyna od dużej wartości i **maleje** (np. c(i+1) = α·c(i), α < 1)
- Rozwiązanie zbliża się do brzegu od wewnątrz
- Tworzy "barierę" na granicy obszaru dopuszczalnego

### Algorytm metody funkcji kary (pseudokod)

```
Dane: x⁽⁰⁾, c⁽¹⁾ > 0, α (skala), ε, Nmax

1. i = 0
2. Powtarzaj:
   i = i + 1
   
   3. Zdefiniuj: F⁽ⁱ⁾(x) = f(x) + c⁽ⁱ⁾·S(x)
   
   4. Znajdź x⁽ⁱ⁾ minimalizując F⁽ⁱ⁾ startując z x⁽ⁱ⁻¹⁾
      (używając np. Nelder-Mead)
   
   5. c⁽ⁱ⁺¹⁾ = α·c⁽ⁱ⁾
   
   6. Jeśli fcalls > Nmax: return error

3. Aż ||x⁽ⁱ⁾ - x⁽ⁱ⁻¹⁾|| < ε
4. Zwróć x* = x⁽ⁱ⁾
```

### Szablon implementacji funkcji kary (`opt_alg.cpp`):

```cpp
solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, 
             double c, double dc, double epsilon, int Nmax, 
             matrix ud1, matrix ud2)
{
    try {
        solution Xopt;
        // Tu wpisz kod funkcji
        // ff - funkcja celu z karą
        // x0 - punkt startowy
        // c - początkowy współczynnik kary
        // dc - współczynnik skalowania (alpha)
        // epsilon - dokładność
        // Nmax - maksymalna liczba wywołań

        return Xopt;
    }
    catch (string ex_info)
    {
        throw ("solution pen(...):\n" + ex_info);
    }
}
```

---

## Problem testowy vs rzeczywisty

### Różnice w podejściu

| Aspekt | Problem testowy | Problem rzeczywisty |
|--------|-----------------|---------------------|
| Funkcja celu | Analityczna (ff3T) | Wynik symulacji (ff3R) |
| Metoda kary | Zewnętrzna i wewnętrzna | Tylko zewnętrzna |
| Punkt startowy (wewn.) | Musi być w obszarze | N/A |
| Liczba optymalizacji | 100 (dla każdego a) | 1 |
| Ograniczenia | 3 nierówności | 1 nierówność (przejście przez punkt) |

### Dlaczego dla problemu rzeczywistego tylko zewnętrzna kara?

1. **Trudność znalezienia punktu startowego** - dla wewnętrznej funkcji kary punkt musi być w obszarze dopuszczalnym, co dla problemu z symulacją może być trudne do zagwarantowania.

2. **Charakter ograniczenia** - ograniczenie "piłka musi minąć punkt (5,50) w odległości ≤ 2m" jest trudne do sprawdzenia bez przeprowadzenia symulacji.

---

## Analiza istniejącej implementacji

### Aktualny stan `lab3()`:

```cpp
void lab3()
{
    double alpha = 1;      // Współczynnik odbicia
    double beta = 0.5;     // Współczynnik zawężenia
    double gamma = 2;      // Współczynnik ekspansji
    double delta = 0.5;    // Współczynnik redukcji
    double epsilon = 1e-3; // Dokładność
    int Nmax = 10000;      // Max wywołań funkcji celu

    int n = 100;           // Liczba optymalizacji (nieużywana)

    solution optX;

    double s = 0.01;       // Bok sympleksu
    matrix k0 = rand_mat(2);

    // Losowanie punktu startowego w [4, 24]
    for (int j = 0; j < 2; ++j)
    {
        k0(j) = 20 * k0(j) + 4;
    }
    
    std::cout << k0(0) << " y: " << k0(1);
    optX = sym_NM(ff3T, k0, s, alpha, beta, gamma, delta, epsilon, Nmax);
    std::cout << "X0: " << k0 << " optimized: " << optX.y;
}
```

### Co należy jeszcze zaimplementować:

1. **Funkcja `pen()`** - metoda funkcji kary
2. **Funkcje kary** (zewnętrzna i wewnętrzna) dla testowej funkcji celu
3. **Funkcja celu dla problemu rzeczywistego** (ff3R) - spadająca piłka
4. **Pętla 100 optymalizacji** z zapisem wyników
5. **Generowanie punktu startowego w obszarze dopuszczalnym** (dla wewnętrznej kary)

### Schemat kompletnej implementacji:

```cpp
// Zewnętrzna funkcja kary dla problemu testowego
matrix ff3T_zewn(matrix x, matrix ud1, matrix ud2) {
    double x1 = m2d(x(0));
    double x2 = m2d(x(1));
    double a = m2d(ud1);    // Parametr ograniczenia (4, 4.4934, lub 5)
    double c = m2d(ud2);    // Współczynnik kary
    
    // Ograniczenia
    double g1 = -x1 + 1;
    double g2 = -x2 + 1;
    double g3 = sqrt(x1*x1 + x2*x2) - a;
    
    // Funkcja kary zewnętrznej
    double S = pow(max(0.0, g1), 2) + 
               pow(max(0.0, g2), 2) + 
               pow(max(0.0, g3), 2);
    
    // Funkcja celu
    double common = M_PI * sqrt(pow(x1/M_PI, 2) + pow(x2/M_PI, 2));
    double f = sin(common) / common;
    
    return matrix(f + c * S);
}

// Wewnętrzna funkcja kary dla problemu testowego
matrix ff3T_wewn(matrix x, matrix ud1, matrix ud2) {
    double x1 = m2d(x(0));
    double x2 = m2d(x(1));
    double a = m2d(ud1);
    double c = m2d(ud2);
    
    double g1 = -x1 + 1;
    double g2 = -x2 + 1;
    double g3 = sqrt(x1*x1 + x2*x2) - a;
    
    // Sprawdzenie czy punkt jest wewnątrz
    if (g1 >= 0 || g2 >= 0 || g3 >= 0) {
        return matrix(1e10);  // Bardzo duża wartość
    }
    
    // Funkcja kary wewnętrznej
    double S = -1.0/g1 - 1.0/g2 - 1.0/g3;
    
    // Funkcja celu
    double common = M_PI * sqrt(pow(x1/M_PI, 2) + pow(x2/M_PI, 2));
    double f = sin(common) / common;
    
    return matrix(f + c * S);
}
```

---

## Podsumowanie algorytmów kary

### Zewnętrzna funkcja kary:
```
ZALETY:
+ Punkt startowy może być gdziekolwiek
+ Łatwa implementacja
+ Dobre dla problemów z trudnym obszarem dopuszczalnym

WADY:
- Rozwiązanie może "oscylować" wokół granicy
- Wymaga dużych wartości c dla dokładnego rozwiązania
```

### Wewnętrzna funkcja kary:
```
ZALETY:
+ Rozwiązanie zawsze w obszarze dopuszczalnym
+ Lepsze zachowanie numeryczne blisko optimum

WADY:
- Wymaga punktu startowego w obszarze
- Problemy numeryczne gdy gᵢ → 0
- Trudniejsza implementacja
```

### Kiedy używać której metody:

| Sytuacja | Zalecana metoda |
|----------|-----------------|
| Łatwy obszar dopuszczalny | Wewnętrzna |
| Trudny obszar dopuszczalny | Zewnętrzna |
| Problem z symulacją | Zewnętrzna |
| Wymagana gwarancja dopuszczalności | Wewnętrzna |
| Szybka konwergencja | Zewnętrzna (z dużym c) |

---

*Dokument przygotowany jako wyjaśnienie kodu do ćwiczenia K3 z przedmiotu Metody Optymalizacji.*
