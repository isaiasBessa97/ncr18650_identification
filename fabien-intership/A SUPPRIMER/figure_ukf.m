% Figure 2 : Transformée Unscented (UKF)
figure('Color', 'w'); hold on;

% Définition de l'axe X et de la distribution d'entrée (Gaussienne bleue)
x = linspace(-12, 12, 300);
pdf_x = normpdf(x, 0, 3) * 30; % Mise à l'échelle pour l'affichage
plot(x, pdf_x, 'b-', 'LineWidth', 2);

% Zone ombrée sous la courbe d'entrée (x_mean +/- sigma)
x_shade = linspace(-4, 4, 100);
y_shade = normpdf(x_shade, 0, 3) * 30;
fill([x_shade, fliplr(x_shade)], [y_shade, zeros(1, 100)], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Fonction non linéaire arbitraire (Courbe verte)
f_nl = @(x) 8 * exp(-(x/6).^2) + 0.15*x + 2;
plot(x, f_nl(x), 'Color', [0 0.5 0], 'LineWidth', 2);

% Les 3 points Sigma traversant la fonction
pts_x = [-6, 0, 7.5];
plot(pts_x, f_nl(pts_x), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);

% Masquer les axes standards pour un aspect purement illustratif
axis off;
text(0, -1, '$\hat{x}_{n+1,n}$', 'Interpreter', 'latex', 'HorizontalAlignment', 'center', 'FontSize', 14);
text(-4.5, -1, '$\hat{x}_{n+1,n} - p_{n+1,n}$', 'Interpreter', 'latex', 'HorizontalAlignment', 'center', 'FontSize', 12);
text(4.5, -1, '$\hat{x}_{n+1,n} + p_{n+1,n}$', 'Interpreter', 'latex', 'HorizontalAlignment', 'center', 'FontSize', 12);