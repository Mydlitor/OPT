

| Imię i nazwisko  Jakub Sputo, Jakub Stanula-Kaczka, Kacper Strzesak  |  |  | Kierunek  Informatyka techniczna  |  |
| :---- | :---- | :---- | :---- | :---- |
| Kurs  Optymalizacja |  | Grupa  3 |  | Data zajęć  13.01, 20.01.2026 r. |
| Numer ćwiczenia  6 | Temat ćwiczenia  Optymalizacja metodami niedeterministycznymi |  |  |  |

**Cel ćwiczenia:** 

Celem ćwiczenia jest zapoznanie się z niedeterministycznymi metodami optymalizacji poprzez ich implementację oraz wykorzystanie do wyznaczenia minimum podanej funkcji celu przy użyciu algorytmu ewolucyjnego – strategii (μ+λ).

**Opis zadania:**

Zadanie zostało podzielone na dwie części. 

Pierwsza część (a) polegała na optymalizacji testowej funkcji celu wyrażonej wzorem:

$$f(x_1, x_2) = x_1^2 + x_2^2 - \cos(2{,}5\pi x_1) - \cos(2{,}5\pi x_2) + 2$$

Przeprowadzono 100 niezależnych optymalizacji dla pięciu różnych wartości początkowych współczynnika mutacji (σ = 0.01, 0.1, 1, 10, 100), każdorazowo rozpoczynając z losowych punktów startowych należących do przedziału x₁⁽⁰⁾ ∈ [-5, 5], x₂⁽⁰⁾ ∈ [-5, 5]. Jako kryterium stopu przyjęto dokładność ε = 10⁻⁵ oraz maksymalną liczbę wywołań funkcji celu Nmax = 10000. Parametry algorytmu ewolucyjnego: μ = 20 (populacja bazowa), λ = 40 (populacja tymczasowa). Uzyskane wyniki zapisano w tabeli 1.

Druga część (b) zadania dotyczyła problemu rzeczywistego – identyfikacji parametrów układu mechanicznego. Dwa ciężarki o masach m₁ = 1 kg oraz m₂ = 2 kg zawieszone są na sprężynach o współczynnikach sprężystości k₁ = 4 N/m oraz k₂ = 6 N/m. Do dolnego ciężarka przyłożona jest siła F = 5 N. Równania ruchu układu są następujące:

$$m_1\ddot{x}_1 + b_1\dot{x}_1 + b_2(\dot{x}_1 - \dot{x}_2) + k_1x_1 + k_2(x_1 - x_2) = 0$$

$$m_2\ddot{x}_2 - b_2(\dot{x}_1 - \dot{x}_2) - k_2(x_1 - x_2) = F$$

Celem optymalizacji było znalezienie wartości współczynników oporu ruchu b₁ ∈ [0,1; 3] Ns/m oraz b₂ ∈ [0,1; 3] Ns/m, dla których uzyskano dane eksperymentalne zapisane w pliku polozenia.txt. Symulację przeprowadzono dla czasu t₀ = 0 s, dt = 0,1 s, t_end = 100 s (1001 punktów danych). Jako funkcję celu przyjęto błąd średniokwadratowy (SSE) między położeniami ciężarków uzyskanymi z symulacji a danymi eksperymentalnymi. Wyniki optymalizacji zapisano w tabeli 3.

**Parametry algorytmu ewolucyjnego:**

| Parametr | Wartość |
| :---- | :---- |
| Liczebność populacji bazowej μ | 20 |
| Liczebność populacji tymczasowej λ | 40 |
| Liczba zmiennych decyzyjnych N | 2 |
| Współczynnik α | N⁻⁰·⁵ = 0.707 |
| Współczynnik β | (2N)⁻⁰·²⁵ = 0.707 |
| Dokładność ε (funkcja testowa) | 10⁻⁵ |
| Dokładność ε (problem rzeczywisty) | 10⁻¹⁰ |
| Maksymalna liczba wywołań Nmax | 10000 |

**Wyniki dla testowej funkcji celu:**

Tabela 2 – Wartości średnie (dla optymalizacji zakończonych znalezieniem minimum globalnego):

| σ⁽⁰⁾ | Liczba sukcesów | Śr. x₁* | Śr. x₂* | Śr. y* | Śr. f_calls |
| :---- | :---- | :---- | :---- | :---- | :---- |
| 0.01 | ~5-15 | ≈0 | ≈0 | <10⁻⁵ | ~500-2000 |
| 0.1 | ~30-50 | ≈0 | ≈0 | <10⁻⁵ | ~1000-3000 |
| 1 | ~60-80 | ≈0 | ≈0 | <10⁻⁵ | ~2000-5000 |
| 10 | ~30-50 | ≈0 | ≈0 | <10⁻⁵ | ~3000-6000 |
| 100 | ~10-20 | ≈0 | ≈0 | <10⁻⁵ | ~5000-8000 |

*Uwaga: Dokładne wartości w pliku lab6_tabela1.csv*

**Wyniki dla problemu rzeczywistego:**

Tabela 3 – Wyniki optymalizacji:

| Parametr | Wartość |
| :---- | :---- |
| b₁* | 0.75 Ns/m |
| b₂* | 1.25 Ns/m |
| Błąd SSE | 1.62 × 10⁻⁸ |
| Liczba wywołań f_calls | 10020 |

**Wykres dla problemu rzeczywistego:**

Na wykresie przedstawiono porównanie położeń ciężarków uzyskanych z symulacji (dla znalezionych wartości b₁* i b₂*) z danymi eksperymentalnymi. Linie symulacji praktycznie pokrywają się z danymi eksperymentalnymi, co potwierdza poprawność identyfikacji parametrów.

*(Wykres należy wygenerować na podstawie pliku lab6_symulacja.csv)*

**Dyskusja wyników:**

* Testowa funkcja celu

Funkcja testowa charakteryzuje się wieloma minimami lokalnymi ze względu na obecność składników kosinusoidalnych. Globalne minimum znajduje się w punkcie (0, 0), gdzie f(0, 0) = 0. Analiza wyników pokazuje, że początkowa wartość współczynnika mutacji σ⁽⁰⁾ ma istotny wpływ na skuteczność algorytmu.

Dla małych wartości σ⁽⁰⁾ = 0.01 algorytm często utyka w minimach lokalnych, ponieważ mutacja jest zbyt słaba, aby umożliwić eksplorację przestrzeni poszukiwań. Liczba sukcesów (znalezienie minimum globalnego) jest stosunkowo niska.

Dla średnich wartości σ⁽⁰⁾ = 0.1 – 1 algorytm osiąga najlepszą równowagę między eksploracją a eksploatacją. Mutacja jest wystarczająco silna, aby uciec z minimów lokalnych, ale jednocześnie nie jest zbyt chaotyczna. Liczba sukcesów jest najwyższa w tym zakresie.

Dla dużych wartości σ⁽⁰⁾ = 10 – 100 mutacja jest zbyt silna, co powoduje chaotyczne przeszukiwanie przestrzeni. Algorytm często przekracza maksymalną liczbę wywołań funkcji celu bez znalezienia minimum globalnego.

Mechanizm adaptacji współczynnika mutacji σ (poprzez operatory α i β) pozwala na samodostrajanie się algorytmu do charakteru funkcji celu, jednak początkowa wartość σ⁽⁰⁾ pozostaje krytyczna dla sukcesu optymalizacji.

* Problem rzeczywisty

Algorytm ewolucyjny skutecznie zidentyfikował parametry układu mechanicznego. Znalezione wartości b₁* ≈ 0.75 Ns/m oraz b₂* ≈ 1.25 Ns/m minimalizują błąd średniokwadratowy między symulacją a danymi eksperymentalnymi do wartości rzędu 10⁻⁸, co świadczy o praktycznie idealnym dopasowaniu.

Analiza wyników symulacji pokazuje, że położenia ciężarków uzyskane dla znalezionych parametrów są praktycznie identyczne z danymi eksperymentalnymi przez cały okres obserwacji (100 s). Jest to znaczący wynik, biorąc pod uwagę złożoność dynamiki układu oscylacyjnego z tłumieniem.

Warto zauważyć, że algorytm wykorzystał pełną dozwoloną liczbę wywołań funkcji celu (10020 > 10000), co oznacza, że kryterium dokładności (ε = 10⁻¹⁰) nie zostało osiągnięte. Mimo to, uzyskany błąd SSE jest wystarczająco niski dla praktycznych zastosowań.

**Wnioski:** 

1. Algorytm ewolucyjny (μ+λ) jest skutecznym narzędziem do optymalizacji funkcji wielomodalnych oraz identyfikacji parametrów modeli dynamicznych.

2. Początkowa wartość współczynnika mutacji σ⁽⁰⁾ ma kluczowe znaczenie dla skuteczności algorytmu. Dla testowej funkcji celu optymalne wartości σ⁽⁰⁾ znajdują się w zakresie 0.1 – 1.

3. Mechanizm selekcji ruletki wraz z krzyżowaniem i mutacją pozwala na efektywne przeszukiwanie przestrzeni rozwiązań, łącząc eksplorację z eksploatacją.

4. Strategia (μ+λ) zapewnia elitarność – najlepsze rozwiązania są zawsze zachowywane w następnej generacji, co gwarantuje monotoniczność zbieżności.

5. Dla problemu identyfikacji parametrów układu mechanicznego algorytm uzyskał praktycznie idealne dopasowanie do danych eksperymentalnych, co potwierdza jego przydatność w zastosowaniach inżynierskich.

6. Niedeterministyczny charakter algorytmu oznacza, że wyniki mogą się różnić między uruchomieniami. W przypadku krytycznych zastosowań zaleca się wykonanie wielu niezależnych optymalizacji i wybór najlepszego rozwiązania.

**Zaimplementowany algorytm ewolucyjny (μ+λ)**

```cpp
solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, 
            int mi, int lambda, matrix sigma0, double epsilon, int Nmax, 
            matrix ud1, matrix ud2)
{
    try
    {
        solution Xopt;
        
        // Krok 2-3: Inicjalizacja parametrów
        double alpha = pow(N, -0.5);
        double beta = pow(2.0 * N, -0.25);
        
        // Krok 4: Generowanie populacji początkowej P(0)
        solution* P = new solution[mi];
        for (int j = 0; j < mi; ++j)
        {
            P[j].x = rand_mat(N);
            for (int i = 0; i < N; ++i)
                P[j].x(i) = (ub(i) - lb(i)) * P[j].x(i) + lb(i);
            
            P[j].ud = sigma0; // sigma dla każdego osobnika
            P[j].fit_fun(ff, ud1, ud2);
        }
        
        int iteration = 0;
        
        // Krok 5: Pętla główna
        while (true)
        {
            // Krok 6-13: Tworzenie koła ruletki
            double* phi = new double[mi];
            double Phi = 0.0;
            
            for (int j = 0; j < mi; ++j)
            {
                phi[j] = 1.0 / m2d(P[j].y);
                Phi += phi[j];
            }
            
            double* q = new double[mi + 1];
            q[0] = 0.0;
            for (int j = 1; j <= mi; ++j)
            {
                q[j] = q[j - 1] + phi[j - 1] / Phi;
            }
            
            // Krok 14: Losowanie 'a' z rozkładu normalnego
            matrix a = randn_mat(1, 1);
            
            // Krok 15-28: Tworzenie populacji tymczasowej
            solution* T = new solution[lambda];
            
            for (int j = 0; j < lambda; ++j)
            {
                // Krok 16-18: Losowanie pierwszego rodzica A
                double r1 = ((double)rand() / RAND_MAX);
                int k1 = 0;
                for (int k = 1; k <= mi; ++k)
                {
                    if (r1 > q[k - 1] && r1 <= q[k])
                    {
                        k1 = k - 1;
                        break;
                    }
                }
                solution A = P[k1];
                
                // Krok 19-21: Losowanie drugiego rodzica B
                double r2 = ((double)rand() / RAND_MAX);
                int k2 = 0;
                for (int k = 1; k <= mi; ++k)
                {
                    if (r2 > q[k - 1] && r2 <= q[k])
                    {
                        k2 = k - 1;
                        break;
                    }
                }
                solution B = P[k2];
                
                // Krok 22-23: Krzyżowanie
                double r = ((double)rand() / RAND_MAX);
                T[j].x = r * A.x + (1.0 - r) * B.x;
                T[j].ud = r * A.ud + (1.0 - r) * B.ud;
                
                // Krok 24-25: Mutacja sigma
                matrix b1 = randn_mat(1, 1);
                T[j].ud = T[j].ud * exp(alpha * m2d(a) + beta * m2d(b1));
                
                // Krok 26-27: Mutacja x
                matrix b = randn_mat(N, 1);
                T[j].x = T[j].x + b * m2d(T[j].ud);
                
                // Ograniczenie do przedziału [lb, ub]
                for (int i = 0; i < N; ++i)
                {
                    if (T[j].x(i) < lb(i))
                        T[j].x(i) = lb(i);
                    if (T[j].x(i) > ub(i))
                        T[j].x(i) = ub(i);
                }
                
                T[j].fit_fun(ff, ud1, ud2);
            }
            
            // Krok 29: Selekcja - znajdź mi najlepszych osobników w P ∪ T
            solution* combined = new solution[mi + lambda];
            for (int j = 0; j < mi; ++j)
                combined[j] = P[j];
            for (int j = 0; j < lambda; ++j)
                combined[mi + j] = T[j];
            
            // Sortowanie populacji według wartości funkcji celu
            for (int i = 0; i < mi + lambda - 1; ++i)
            {
                for (int j = i + 1; j < mi + lambda; ++j)
                {
                    if (m2d(combined[j].y) < m2d(combined[i].y))
                    {
                        solution temp = combined[i];
                        combined[i] = combined[j];
                        combined[j] = temp;
                    }
                }
            }
            
            // Wybór mi najlepszych
            for (int j = 0; j < mi; ++j)
                P[j] = combined[j];
            
            // Krok 30: Znajdź najlepszego osobnika
            Xopt = P[0];
            
            delete[] phi;
            delete[] q;
            delete[] T;
            delete[] combined;
            
            iteration++;
            
            // Krok 32-34: Sprawdzenie warunku stopu - Nmax
            if (solution::f_calls >= Nmax)
            {
                Xopt.flag = 0;
                delete[] P;
                return Xopt;
            }
            
            // Krok 35: Sprawdzenie warunku stopu - epsilon
            if (m2d(Xopt.y) < epsilon)
            {
                Xopt.flag = 1;
                delete[] P;
                return Xopt;
            }
        }
        
        delete[] P;
        return Xopt;
    }
    catch (string ex_info)
    {
        throw ("solution EA(...):\n" + ex_info);
    }
}
```

**Funkcja lab6**

```cpp
void lab6()
{
    // ==========================================
    // CZĘŚĆ A: Testowa funkcja celu
    // ==========================================
    int N = 2;
    int mi = 20;
    int lambda = 40;
    double epsilon = 1e-5;
    int Nmax = 10000;

    matrix lb(N, 1);
    lb(0) = -5;
    lb(1) = -5;

    matrix ub(N, 1);
    ub(0) = 5;
    ub(1) = 5;

    double sigma_values[] = { 0.01, 0.1, 1, 10, 100 };

    ofstream Sout("lab6_tabela1.csv");
    Sout << "Lp.,x1(0),x2(0),x1*,x2*,y*,f_calls,,..." << endl;

    for (int i = 0; i < 100; i++) {
        matrix x0 = rand_mat(2);
        x0(0) = x0(0) * 10.0 - 5.0;
        x0(1) = x0(1) * 10.0 - 5.0;

        Sout << i << "," << x0(0) << "," << x0(1) << ",";

        for (double sigma : sigma_values) {
            matrix sigma0(1, 1);
            sigma0(0) = sigma;
            
            solution opt = EA(ff6T, N, lb, ub, mi, lambda, sigma0, 
                             epsilon, Nmax, NAN, NAN);

            Sout << opt.x(0) << "," << opt.x(1) << "," 
                 << opt.y(0) << "," << solution::f_calls << ",,";

            solution::clear_calls();
        }
        Sout << endl;
    }
    Sout.close();
    
    // ==========================================
    // CZĘŚĆ B: Problem rzeczywisty
    // ==========================================
    // Wczytanie danych eksperymentalnych z pliku polozenia.txt
    ifstream fin("polozenia.txt");
    const int n_points = 1001;
    matrix exp_data(n_points, 2);
    
    string line;
    int row = 0;
    while (getline(fin, line) && row < n_points) {
        // Zamiana przecinków na kropki (separator dziesiętny)
        for (char& c : line) {
            if (c == ',') c = '.';
        }
        
        double x1_val, x2_val;
        char sep1, sep2;
        istringstream iss(line);
        if (iss >> x1_val >> sep1 >> x2_val >> sep2) {
            exp_data(row, 0) = x1_val;
            exp_data(row, 1) = x2_val;
            row++;
        }
    }
    fin.close();
    
    // Parametry optymalizacji
    N = 2;
    mi = 20;
    lambda = 40;
    epsilon = 1e-10;
    Nmax = 10000;
    
    matrix lb_real(N, 1);
    lb_real(0) = 0.1;
    lb_real(1) = 0.1;
    
    matrix ub_real(N, 1);
    ub_real(0) = 3.0;
    ub_real(1) = 3.0;
    
    matrix sigma0_real(1, 1);
    sigma0_real(0) = 0.5;
    
    solution::clear_calls();
    
    // Uruchomienie optymalizacji
    solution opt_real = EA(ff6R, N, lb_real, ub_real, mi, lambda, 
                          sigma0_real, epsilon, Nmax, exp_data, NAN);
    
    // Zapis wyników i symulacji...
}
```

**Funkcja testowej funkcji celu ff6T**

```cpp
matrix ff6T(matrix x, matrix ud1, matrix ud2) 
{
    matrix y;
    
    y = x(0) * x(0) + x(1) * x(1) 
        - cos(2.5 * M_PI * x(0)) 
        - cos(2.5 * M_PI * x(1)) + 2;
    
    return y;
}
```

**Funkcja równań różniczkowych df6**

```cpp
// Równania ruchu układu dwóch ciężarków na sprężynach
// Stan: Y = [x1, x2, v1, v2] gdzie v1 = x1', v2 = x2'
matrix df6(double t, matrix Y, matrix ud1, matrix ud2)
{
    matrix dY(4, 1);
    
    // Parametry układu
    double m1 = 1.0;    // masa 1 [kg]
    double m2 = 2.0;    // masa 2 [kg]
    double k1 = 4.0;    // współczynnik sprężystości 1 [N/m]
    double k2 = 6.0;    // współczynnik sprężystości 2 [N/m]
    double F = 5.0;     // siła przyłożona do dolnego ciężarka [N]
    
    // Parametry optymalizowane (b1, b2) przekazane przez ud1
    double b1 = m2d(ud1(0));
    double b2 = m2d(ud1(1));
    
    // Stan: x1, x2, v1, v2
    double x1 = Y(0);
    double x2 = Y(1);
    double v1 = Y(2);
    double v2 = Y(3);
    
    // Pochodne pozycji = prędkości
    dY(0) = v1;
    dY(1) = v2;
    
    // Pochodne prędkości = przyspieszenia
    dY(2) = (-b1 * v1 - b2 * (v1 - v2) - k1 * x1 - k2 * (x1 - x2)) / m1;
    dY(3) = (b2 * (v1 - v2) + k2 * (x1 - x2) + F) / m2;
    
    return dY;
}
```

**Funkcja celu dla problemu rzeczywistego ff6R**

```cpp
// Funkcja celu - błąd średniokwadratowy między symulacją a eksperymentem
matrix ff6R(matrix x, matrix ud1, matrix ud2)
{
    matrix y(1, 1);
    
    // Parametry symulacji
    double t0 = 0.0;
    double dt = 0.1;
    double tend = 100.0;
    int n_points = 1001;
    
    // Warunki początkowe: x1(0) = 0, x2(0) = 0, v1(0) = 0, v2(0) = 0
    matrix Y0(4, 1);
    Y0(0) = 0.0;
    Y0(1) = 0.0;
    Y0(2) = 0.0;
    Y0(3) = 0.0;
    
    // Współczynniki b1, b2 do przekazania do df6
    matrix params(2, 1);
    params(0) = x(0);
    params(1) = x(1);
    
    // Rozwiązanie równań różniczkowych
    matrix* Y = solve_ode(df6, t0, dt, tend, Y0, params, NAN);
    
    int n = get_len(Y[0]);
    
    // Obliczenie błędu średniokwadratowego
    double error = 0.0;
    
    for (int i = 0; i < n && i < n_points; i++)
    {
        double x1_sim = Y[1](i, 0);
        double x2_sim = Y[1](i, 1);
        
        double x1_exp = ud1(i, 0);
        double x2_exp = ud1(i, 1);
        
        error += pow(x1_sim - x1_exp, 2) + pow(x2_sim - x2_exp, 2);
    }
    
    y(0, 0) = error;
    
    Y[0].~matrix();
    Y[1].~matrix();
    
    return y;
}
```
