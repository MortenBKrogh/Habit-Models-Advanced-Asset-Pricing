function [percentage] = percentage(vec)

percentage = sum(vec(:)==1) / length(vec);

end