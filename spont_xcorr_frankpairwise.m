% % % %  version 5, try subtract baseline for center of mass

%%

clear all
bin_size=0.01;
xc_len=round(2.00/bin_size);

artifact_rate=100*bin_size;
a=1;
% nr_130_00 = 300;
% nr_131_00 = 660;

% filenames091026
% filenames091112
% filenames091221
% filenames100118
% filenames100125
% filenames100128
% filenames100224
% filenames1003031
% filenames100306
% filenames100319

% filenames100320
% filenames100321
% filenames100324
% filenames100407
% filenames100323
% filenames100325
% filenames100409
% filenames100410
% filenames100421
% filenames100428
% filenames100326
% filenames100506
% filenames1005111
% filenames100514
% filenames1005112
% filenamestotal_4to5
% filenames100519_old_total
% filenames1005TOTAL
% modified_new
filenames100602
% filenames100613
% filenames100619
% filenames100621
% filenames100625
% 
% filenames100629
% filenames1006030

% filenames100702
% filenames100708
xc_pre_all=[];
xc_post_all=[];
xc_cond_all = [];
xc_post1_all=[];
xc_post2_all=[];
xc_post3_all=[];
distance_all=[];
n_xc=0;

%% time
time=datestr(now,30);
fid = fopen([ time '.txt'], 'wt');
fprintf(fid, '%s\n', '*********Selected file**********');

for kk=1:a-1;

load(filenames(kk).name) 
fid = fopen([ time '.txt'], 'at');
fprintf(fid, '%s\n', filenames(kk).name);
    template_st=round(0/bin_size);
    % template_end=round(1.2/bin_size);
    template_baseline_end=round(3/bin_size);

    % template_st=round(0.01/bin_size);
    % template_end=round(1.5/bin_size);


    clear tmp tmp1 tmp2 tmp3 data_pre data_post data_cond
    Datatype = input('Is your data spike(s) or LFP(l) ? ...', 's');
     if Datatype =='s' 
  %% Spike 
  for ii=[1:16]
        tmp=num2str(ii);
        if length(tmp)==1
            tmp=['0' tmp];
        end
        var_name=['nr_0' tmp '_00'];
        tmp=eval(var_name);
 
        tmp1=hist(tmp, [0:bin_size:min(nr_130_00)]);
        data_pre(ii,:)=tmp1(1:end-1);

         tmp2=hist(tmp,[max(nr_131_00):bin_size:900]);
        data_post(ii,:)=tmp2(2:end-1);
    
        tmp3=hist(tmp, [min(nr_130_00):bin_size:max(nr_131_00)+wave_dur]);
         data_cond(ii,:)=tmp3(2:end-1);
    
  end
    data_pre_raw=data_pre(included_channel(kk).channel,:);
    data_post_raw=data_post(included_channel(kk).channel,:);
    data_cond_raw=data_cond(included_channel(kk).channel,:);

    data_pre=data_pre_raw; art_ind=find(mean(data_pre,1)>artifact_rate | prod(data_pre,1)>0); data_pre(:,art_ind)=0;
%     data_pre=data_pre_raw; art_ind=find(sum(data_pre,1)>artifact_rate); data_pre(:,art_ind)=0;
    data_cond=data_cond_raw; art_ind=find(mean(data_cond,1)>artifact_rate | prod(data_cond,1)>0); data_cond(:,art_ind)=0;
    data_post=data_post_raw; art_ind= find(mean(data_post,1)>artifact_rate | prod(data_post,1)>0); data_post(:,art_ind)=0;
    data_post1=data_post(:,1:round(20/bin_size));
    data_post2=data_post(:,round(20/bin_size)+1:round(40/bin_size));
    data_post3=data_post(:,round(40/bin_size)+1:min(round((60/bin_size)), size(data_post,2)));

 elseif Datatype =='l' 

 %% LFP 
 nr_130_00 = 300;nr_131_00 = 660;
  for ii=[1:16]
        tmp=num2str(ii);
        var_name = ['ad_' tmp '_1_kS_s'];
        tmp = eval(var_name);
        tmp1 = tmp(1:300000);
        data_pre(ii,:)=tmp1(1:end-1);
        tmp2= tmp(660000:960000);
        data_post(ii,:)=tmp1(2:end-1);
        tmp3= tmp(300000:660000); 
        data_cond(ii,:)=tmp3(2:end-1);
   end
    data_pre_raw=data_pre(included_channel(kk).channel,:);
    data_post_raw=data_post(included_channel(kk).channel,:);
    data_cond_raw=data_cond(included_channel(kk).channel,:);
    
    data_pre=data_pre_raw;
    data_cond=data_cond_raw;
    data_post=data_post_raw;
    else disp ('Error. Sorry what is your data type?') 
         break 
    end
 
    n_channel=size(data_pre_raw,1);
    RF_position=[1:n_channel]/n_channel;

 %%%%%%%%%%%%%%%%% measure direction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear tmp tmp2 event_raster_cond_stim 

    for ii=1:length(nr_130_00)
          tmp=floor((nr_130_00(ii)-nr_130_00(1))/bin_size);
          if tmp+template_baseline_end<=size(data_cond,2)
            event_raster_cond_stim(ii,:,:)=data_cond(:,tmp+template_st+1:min(tmp+template_baseline_end, size(data_cond,2)));
          end
    end
    cond_template_long=squeeze(sum(event_raster_cond_stim,1));
%% ·­×ª
% cond_template_long=flipud(cond_template_long);

%%%%%%%%%%% subtract baseline %%%%%%%%%%%%%%%%%%%%
    cond_template_long=cond_template_long-repmat(mean(cond_template_long(:, template_end(kk)-template_st+1:end),2),1,size(cond_template_long,2));
    cond_template_long=max(cond_template_long,0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);clf

imagesc(cond_template_long)
title(['cond-template-long' num2str(kk)])
% set(gca,'YTick',[]);
xlabel('time(ms)');ylabel('channel');
saveppt2('602','figure','1');
cond_template=cond_template_long(:,1:template_end(kk));
% 
%%%%%%%%%%%%%%%%%% re-sequence %%%%%%%%%%%%%%%%%%%%%%%%%%%
% [channel_seq_template  sorted_COM_tmplate wave_speed_template wave_CC_template]=wave_dir_ver3(cond_template);% sorting sequrnce
channel_seq_template=[1:size(data_pre,1)]; % sequence without sorting
% 
channel_seq_template
cond_template=cond_template(channel_seq_template,:);
% % 
data_pre(:,:,1)=data_pre;
data_cond(:,:,1)=data_cond;
data_post(:,:,1)=data_post;
% data_post1(:,:,1)=data_post1;
% data_post2(:,:,1)=data_post2;
% data_post3(:,:,1)=data_post3;
% 
% % data_pre(:,:,1)=data_pre(channel_seq_template,:);
% % data_cond(:,:,1)=data_cond(channel_seq_template,:);
% % data_post(:,:,1)=data_post(channel_seq_template,:);
% % data_post1(:,:,1)=data_post1(channel_seq_template,:);
% % data_post2(:,:,1)=data_post2(channel_seq_template,:);
% % data_post3(:,:,1)=data_post3(channel_seq_template,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [tmp  sorted_COM_tmplate_tmp wave_speed_template wave_CC_template]=wave_dir_ver3(cond_template);
figure(2);clf
subplot(2,1,1)
imagesc(cond_template)
title(['cond-template' num2str(kk)])
subplot(2,1,2)
imagesc(cond_template./repmat(max(cond_template,[],2),1,size(cond_template,2)))
title(['cond-template-norm' num2str(kk)])
hold on
% subplot(3,1,2)
% imagesc(flipud(cond_template));
% title(['cond-template_descend' num2str(kk)])
% plot(sorted_COM_tmplate, [1:length(sorted_COM_tmplate)],'ro')
saveppt2('602','figure','2');
% 

    t=[-xc_len:xc_len];
% % % %     for channel_shift=channel_shift_rg
% % % %         xc_pre(kk,channel_shift,:)=xc_between_channel(data_pre, channel_shift, xc_len, 0);
% % % %         xc_post(kk,channel_shift,:)=xc_between_channel(data_post, channel_shift, xc_len, 0);
% % % %         xc_post1(kk,channel_shift,:)=xc_between_channel(data_post1, channel_shift, xc_len, 0);
% % % %         xc_post2(kk,channel_shift,:)=xc_between_channel(data_post2, channel_shift, xc_len, 0);
% % % %         xc_post3(kk,channel_shift,:)=xc_between_channel(data_post3, channel_shift, xc_len, 0);
% % % %     end

    [xc_pre distance]=xc_vs_RFdistance(data_pre, RF_position, xc_len);
    [xc_post distance]=xc_vs_RFdistance(data_post, RF_position, xc_len);
    [xc_cond distance]=xc_vs_RFdistance(data_cond, RF_position, xc_len);
%     [xc_post1 distance]=xc_vs_RFdistance(data_post1, RF_position, xc_len);
%     [xc_post2 distance]=xc_vs_RFdistance(data_post2, RF_position, xc_len);
%     [xc_post3 distance]=xc_vs_RFdistance(data_post3, RF_position, xc_len);
    
    xc_pre_all=[xc_pre_all;squeeze(xc_pre)];
    xc_post_all=[xc_post_all;squeeze(xc_post)];
    xc_cond_all=[xc_cond_all;squeeze(xc_cond)];
%     xc_post1_all=[xc_post1_all;squeeze(xc_post1)];
%     xc_post2_all=[xc_post2_all;squeeze(xc_post2)];
%     xc_post3_all=[xc_post3_all;squeeze(xc_post3)];
    distance_all=[distance_all distance];
    
    n_xc=n_xc+size(xc_pre,2);


end
fclose(fid)
[tmp, RF_rank]=sort(distance_all(:));

% figure(3);clf
% subplot(4,1,1)
% imagesc(xc_pre_all(RF_rank,:))
% hold on
% plot([xc_len+1 xc_len+1], [0 n_xc],'k')
% title('pre')
% subplot(4,1,2)
% imagesc(xc_post_all(RF_rank,:))
% hold on
% plot([xc_len+1 xc_len+1], [0 n_xc],'k')
% title('post')
% subplot(4,1,3)
% imagesc(xc_post_all(RF_rank,:)-xc_pre_all(RF_rank,:))
% hold on
% plot([xc_len+1 xc_len+1], [0 n_xc],'k')
% title('diff')
% subplot(4,1,4)
% tmp=mean(xc_pre_all,1); tmp=tmp-min(tmp); tmp=tmp/max(tmp); plot(t,tmp)
% hold on
% tmp=mean(xc_post_all,1); tmp=tmp-min(tmp); tmp=tmp/max(tmp); plot(t,tmp,'r')
% plot([0 0],[0 1],'g:')
% saveppt2('602','figure','3');
% clear xc_pre_m xc_post_m
dist_boundary=[0:10]/10;
nonzero_bin=[];
for ii=1:length(dist_boundary)-1
    ind=find(distance_all>dist_boundary(ii) & distance_all<=dist_boundary(ii+1));
    if length(ind)>0
        xc_pre_m(ii,:)=mean(xc_pre_all(ind,:),1);
        xc_post_m(ii,:)=mean(xc_post_all(ind,:),1);
        xc_cond_m(ii,:)=mean(xc_cond_all(ind,:),1);
        nonzero_bin=[nonzero_bin ii];
    else
        xc_pre_m(ii,:)=zeros(1,2*xc_len+1);
        xc_post_m(ii,:)=zeros(1,2*xc_len+1);
        xc_cond_m(ii,:)=zeros(1,2*xc_len+1);
    end

    tmp=xc_pre_m(ii,:); tmp=tmp-min(tmp); xc_pre_m_COM(ii)=sum(tmp.*[1:2*xc_len+1])/sum(tmp);
    tmp=xc_post_m(ii,:); tmp=tmp-min(tmp); xc_post_m_COM(ii)=sum(tmp.*[1:2*xc_len+1])/sum(tmp);
    tmp=xc_cond_m(ii,:); tmp=tmp-min(tmp); xc_cond_m_COM(ii)=sum(tmp.*[1:2*xc_len+1])/sum(tmp);
    tmp=xc_post_m(ii,:)-xc_pre_m(ii,:); tmp=tmp-min(tmp); xc_diff_m_COM(ii)=sum(tmp.*[1:2*xc_len+1])/sum(tmp);
    
end

%%%%%%%% a bunch of manipulations %%%%%%%%%%%%%%%%%%%%%%%
xc_pre_m(:,max(t)+1)=(xc_pre_m(:,max(t))+xc_pre_m(:,max(t)+2))/2;
xc_post_m(:,max(t)+1)=(xc_post_m(:,max(t))+xc_post_m(:,max(t)+2))/2;
xc_cond_m(:,max(t)+1)=(xc_cond_m(:,max(t))+xc_cond_m(:,max(t)+2))/2;

xc_pre_m=xc_pre_m./repmat(sum(xc_pre_m,2),1,size(xc_pre_m,2));
xc_post_m=xc_post_m./repmat(sum(xc_post_m,2),1,size(xc_post_m,2));
xc_cond_m=xc_cond_m./repmat(sum(xc_cond_m,2),1,size(xc_cond_m,2));
nonzero_bin=nonzero_bin(1:end-1);
%%%%%%%% a bunch of manipulations %%%%%%%%%%%%%%%%%%%%%%%


figure(4);clf

subplot(4,1,1)
imagesc(flipud(xc_pre_m(nonzero_bin,:)))
% xlabel('time(s)')
% ylabel('channel distance')
set(gca,'YTick',[]);set(gca,'XTick',[]);

hold on
plot([xc_len+1 xc_len+1], [0 n_xc],'k')
plot(fliplr(xc_pre_m_COM(nonzero_bin)),nonzero_bin-min(nonzero_bin)+1,'w')
title('pre')
subplot(4,1,2)
imagesc(flipud(xc_post_m(nonzero_bin,:)))
set(gca,'YTick',[]);set(gca,'XTick',[]);
hold on
plot([xc_len+1 xc_len+1], [0 n_xc],'k')
plot(fliplr(xc_post_m_COM(nonzero_bin)),nonzero_bin-min(nonzero_bin)+1,'w')
title('post')

subplot(4,1,3)
imagesc(flipud(xc_cond_m(nonzero_bin,:)))
% xlabel('time(s)')
% ylabel('channel distance')
set(gca,'YTick',[]);set(gca,'XTick',[]);
hold on
plot([xc_len+1 xc_len+1], [0 n_xc],'k')
plot(fliplr(xc_cond_m_COM(nonzero_bin)),nonzero_bin-min(nonzero_bin)+1,'w')
title('cond')

subplot(4,1,4)
imagesc(flipud(xc_post_m(nonzero_bin,:)-xc_pre_m(nonzero_bin,:)))
% imagesc(flipud(xc_C_post_m/mean(xc_C_post_m(~isnan(xc_C_post_m)))-xc_C_pre_m/mean(xc_C_pre_m(~isnan(xc_C_pre_m)))))
% imagesc(flipud(xc_C_post_m/max(xc_C_post_m(:))-xc_C_pre_m/max(xc_C_pre_m(:))))

set(gca,'YTick',[]);set(gca,'XTick',[]);
hold on
plot([xc_len+1 xc_len+1], [0 n_xc],'k')
plot(fliplr(xc_diff_m_COM(nonzero_bin)),nonzero_bin-min(nonzero_bin)+1,'w')
title('diff')
xlabel('time(s)')
% ylabel('channel distance')
saveppt2('602','figure','4');

