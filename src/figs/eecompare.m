load('../solutions/eeerrors_simfp7ART')
mean7ART = R.meanEE;
max7ART = R.maxEE;
load('../solutions/eeerrors_simfp9ART')
mean9ART = R.meanEE;
max9ART = R.maxEE;
load('../solutions/eeerrors_simfp11ART')
mean11ART = R.meanEE;
max11ART = R.maxEE;
load('../solutions/eeerrors_simfp13ART')
mean13ART = R.meanEE;
max13ART = R.maxEE;
load('../solutions/eeerrors_simfp7Gust')
mean7Gust = R.meanEE;
max7Gust = R.maxEE;
load('../solutions/eeerrors_simfp9Gust')
mean9Gust = R.meanEE;
max9Gust = R.maxEE;
load('../solutions/eeerrors_simfp11Gust')
mean11Gust = R.meanEE;
max11Gust = R.maxEE;
load('../solutions/eeerrors_simfp13Gust')
mean13Gust = R.meanEE;
max13Gust = R.maxEE;
load('../solutions/eeerrors_simfp15Gust')
mean15Gust = R.meanEE;
max15Gust = R.maxEE;
load('../solutions/eeerrors_simfp17Gust')
mean17Gust = R.meanEE;
max17Gust = R.maxEE;

disp(['ART Consumption Euler: ', num2str(mean7ART(1)),', ',num2str(mean9ART(1)),', ',num2str(mean11ART(1)),', ',num2str(mean13ART(1))])
disp(['ART Firm Pricing: ', num2str(mean7ART(2)),', ',num2str(mean9ART(2)),', ',num2str(mean11ART(2)),', ',num2str(mean13ART(2))])
disp(['ART differences: ',num2str(mean9ART(1)-mean7ART(1)),', ',num2str(mean11ART(1)-mean9ART(1)),', ',num2str(mean13ART(1)-mean11ART(1))])
disp(['Gust Consumption Euler: ', num2str(mean7Gust(1)),', ',num2str(mean9Gust(1)),', ',num2str(mean11Gust(1)),', ',num2str(mean13Gust(1)),', ',num2str(mean15Gust(1)),', ',num2str(mean17Gust(1))])
disp(['Gust Firm Pricing: ', num2str(mean7Gust(2)),', ',num2str(mean9Gust(2)),', ',num2str(mean11Gust(2)),', ',num2str(mean13Gust(2))])
disp(['Gust differences: ',num2str(mean9Gust(1)-mean7Gust(1)),', ',num2str(mean11Gust(1)-mean9Gust(1)),', ',num2str(mean13Gust(1)-mean11Gust(1))])

% disp('Mean Euler Equation Error')
% disp(R.meanEE)
% disp('Max Euler Equation Error')
% disp(R.maxEE)

