function DisplayImage(figure_nr, x_vec, y_vec, Reconstruction, dyn_range, ImageTitle, FontSize)

figure(figure_nr)
imagesc(x_vec, y_vec, 20*log10(abs(Reconstruction)/max(max(abs(Reconstruction)))),[-dyn_range 0])
axis xy image;
h=title(ImageTitle);
set(h, 'FontSize', FontSize);
h = xlabel('x [m]');
set(h, 'FontSize', FontSize);
h = ylabel('y [m]');
set(h, 'FontSize', FontSize);
colorbar('hide');
colormap ('hot');
set(gca, 'FontSize', FontSize);
