from scipy.constants import physical_constants
import matplotlib.pyplot as plt
import scipy.special as sp
import seaborn as sns
import numpy as np

def radial_function(n: int, l: int, r: np.ndarray, a0: float):

    # n --> główna liczba kwantowa
    # l --> poboczna liczba kwantowa
    # r --> współrzędna promieniowa
    # a0 --> promień atomu Bohr'a

    #Wielomian Laguerre'a
    laguerre = sp.genlaguerre(n-l-1, 2*l+1)
    p = 2*r / (n*a0)
    sqrt_part = np.sqrt(((2 / n * a0) ** 3 * (sp.factorial(n - l - 1))) /(2 * n * (sp.factorial(n + l))))

    return sqrt_part * np.exp(-p / 2) * (p ** l) * laguerre(p)

# print(radial_function(1,0,5,0.3))


def angular_function(l: int, m: int, theta: np.ndarray, phi: int):

    # l --> poboczna liczba kwantowa
    # m --> magnetyczna liczba kwantowa
    # theta --> kąt biegunowy
    # phi --> kąt azymutalny

    #Wielomian Legendre'a
    legendre = sp.lpmv(m, l, np.cos(theta))
    sqrt_part = ((-1) ** m) * np.sqrt(((2 * l + 1) * sp.factorial(l - np.abs(m))) /(4 * np.pi * sp.factorial(l + np.abs(m))))
    return sqrt_part * legendre * np.real(np.exp(1.j * m * phi))

# print(angular_function(0,0,5,5))

def wavefunction(n, l, m, a0_scale_factor):

    # skalowanie promienia atomu Bohr'a
    a0 = a0_scale_factor * physical_constants['Bohr radius'][0] * 1e+12

    grid_extent = 480
    grid_resolution = 680
    # x i z, żeby mieć przekrój
    x = z = np.linspace(-grid_extent, grid_extent, grid_resolution)
    x, z = np.meshgrid(x, z)
    #epsilon zeby uniknac dzielenia przez 0 - epsilon reprezentuje najmniejszą możliwą liczbę zmiennoprzecinkową
    eps = np.finfo(float).eps

    r = np.sqrt((x ** 2 + z ** 2))
    theta = np.arctan(z / (x + eps))
    #phi = 0 i jest stałe, ponieważ robimy przekrój w 2D
    phi = 0


    psi = radial_function(n, l, r, a0) * angular_function(l, m, theta, phi)

    return psi



def probability_density(psi):

    #probability density to kwadrat wartości bezwględnej z wavefunction
    return np.abs(psi) ** 2


def visual(n, l, m, a0_scale_factor, colormap='rocket'):

    try:
        sns.color_palette(colormap)
    except ValueError:
        raise ValueError(f'{colormap} is not a recognized Seaborn colormap.')

    plt.rcParams['font.family'] = 'STIXGeneral'
    plt.rcParams['xtick.major.width'] = 4
    plt.rcParams['ytick.major.width'] = 4
    plt.rcParams['xtick.major.size'] = 15
    plt.rcParams['ytick.major.size'] = 15
    plt.rcParams['xtick.labelsize'] = 30
    plt.rcParams['ytick.labelsize'] = 30
    plt.rcParams['axes.linewidth'] = 4

    fig, ax = plt.subplots(figsize=(16, 16.5))
    plt.subplots_adjust(top=0.82)
    plt.subplots_adjust(right=0.97)
    plt.subplots_adjust(left=-0.5)

    psi = wavefunction(n, l, m, a0_scale_factor)
    prob_density = probability_density(psi)

    im = ax.imshow(np.sqrt(prob_density).T, cmap=sns.color_palette(colormap, as_cmap=True))

    cbar = plt.colorbar(im, fraction=0.046, pad=0.03)
    cbar.set_ticks([])
    background_color = sorted(sns.color_palette(colormap, n_colors=100),key=lambda color: 0.2126 * color[0] + 0.7152 * color[1] + 0.0722 * color[2])[0]
    plt.rcParams['text.color'] = '#dfdfdf'
    title_color = '#dfdfdf'
    fig.patch.set_facecolor(background_color)
    cbar.outline.set_visible(False)
    ax.tick_params(axis='x', colors='#c4c4c4')
    ax.tick_params(axis='y', colors='#c4c4c4')
    for spine in ax.spines.values():
        spine.set_color('#c4c4c4')
    ax.set_title('Wizualizacja Orbitali Atomu Wodoru',
                 pad=130, fontsize=44, loc='left', color=title_color)
    ax.text(30, 615, r'$({0}, {1}, {2})$'.format(n, l, m), color='#dfdfdf', fontsize=42)
    ax.text(705, 700, 'Wyższe\nprawd.', fontsize=24)
    ax.text(705, -60, 'Niższe\npraw.', fontsize=24)
    ax.invert_yaxis()

    plt.show()
    # print(prob_density)
    # print(prob_density.shape)
    # print(np.sqrt(prob_density).T)


if __name__ == "__main__":
    running = True
    print("Uruchomiono projekt \"Wizualizacja Orbitali Atomu Wodoru\"")

    while running:
        print("Proszę poprawnie wpisać wartości liczb kwantowych, aby otrzymać wizualizację orbitali atomu wodoru, bądź wpisać \"exit\", aby wyłączyc program:\n")


        #zabezpieczyc petle



        quantum_check = True
        exit = False
        while quantum_check:
            n = input("n - główna liczba kwantowa: ")
            if n == "exit":
                quantum_check = False
                exit = True
                continue
            else:
                try:
                    n = int(n)
                    if n < 1:
                        print("n nie może być mniejsze niż 1\n")
                        continue
                except ValueError:
                    print("n musi być liczbą naturalną większą od 0\n")
                    continue

            l = input("l - poboczna liczba kwantowa: ")
            if l == "exit":
                quantum_check = False
                exit = True
                continue
            else:
                try:
                    l = int(l)
                    if not (0 <= l < n):
                        print("l nie może być mniejsze niż 0 lub większe niż n\n")
                        continue
                except ValueError:
                    print("l musi być liczbą naturalną większą lub równą 0 oraz mniejszą niż n\n")
                    continue

            m = input("m - magnetyczna liczba kwantowa: ")
            if m == "exit":
                quantum_check = False
                exit = True
                continue
            else:
                try:
                    m = int(m)
                    if not (-l <= m <= l):
                        print("m nie może być mniejsze niż -l lub większe niż l\n")
                        continue
                except ValueError:
                    print("m musi być liczbą naturalną w przedziale [-l, l]\n")
                    continue

            a = input("a0_scaling - skalowanie promienia atomu Bohr'a: ")
            if m == "exit":
                quantum_check = False
                exit = True
                continue
            else:
                try:
                    a = float(a)
                    if a <= 0:
                        print("liczba skalująca nie może być mniejsza lub równa 0\n")
                        continue
                except ValueError:
                    print("a0 musi być dodatnią liczbą różną od 0\n")
                    continue

            quantum_check = False



        if exit:
            print("Dziękujemy za skorzystanie z naszego programu!")
            break

        print()
        print("---------------------------------------------------------")
        print()
        visual(n,l,m,a,"mako")








