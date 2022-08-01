clear;clc

video1= VideoReader('microstructure evolution_twinning.avi');
vid = read(video1);
frames = video1.NumberOfFrames;

video2= VideoWriter('test2.avi'); 
video2.FrameRate = 5; video2.Quality= 100; open(video2);

% for x = 50:70
%     imwrite(vid(:,:,:,x),strcat('frame-',num2str(x),'.tiff'));
% end

for x = 150 : frames

      f=figure('visible','off'); imshow(vid(:,:,:,x));
      frame = getframe(gcf); writeVideo(video2,frame);
      f;clf;
end

for k = 1:70

   f=figure('visible','off'); 
   imshow(vid(:,:,:,k));
   frame = getframe(gcf);
   writeVideo(video2,frame);
   f;clf;
end

close(video2);

