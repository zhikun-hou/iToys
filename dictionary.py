from embedder import Embedder,DATASET_ROOT
import translators as ts
import pandas as pd
import os

# 用于NLP中生成双语词典



# 如果两种语言能够互译，说明是高质量单词，放入双语词典中，用作之后的训练和测试
# Tencent AI Lab的中英文词库都各有两百万个词......没法一个一个翻译
# 必须要拼接起来进行翻译
def generate(source_language,target_language,root=DATASET_ROOT):
    SRC_embedder = Embedder.load(source_language,root)
    
    # 将许多个词用换行符拼接起来，再调用百度翻译的api进行查询
    # 不然一个一个查的话，词库里有几百万个词
    forward_words      = SRC_embedder.encoder.keys()
    forward_translate  = query_words(forward_words,source_language,target_language)
    backward_words     = forward_translate.values()
    backward_translate = query_words(backward_words,source_language,target_language)

    # 检查是否能够回译

    DICT_output = dict()
    for raw in forward_words:
        forward  = forward_translate[raw]
        backward = backward_translate[forward]
        if(backward==raw):
            DICT_output[raw] = forward

    # 输出双语词典

    to_path = os.path.join(root,"dict_{}_to_{}.pkl".format(source_language,target_language))
    pd.to_pickle(DICT_output,to_path)
    
    return DICT_output

# 仅限模块内调用，不应该在外部使用
def query_words(word_list,source_language,target_language):
    query = ""
    all_query = []
    # 构造若干次查询
    for word in word_list:
        if(len(query+word+"\n")<5000): # 字数太多会超出翻译上限
            query += word+"\n"
        else:
            all_query.append(query)
            query = word+"\n"
    # 执行若干次查询
    mapping = dict()
    for query in all_query:
        resposne = ts.translate_text(query,"baidu",from_lang=source_language,to_lang=target_language,if_use_preacceleration=False)
        results = resposne.split("\n") # 使用分隔符得到每个词对应的翻译结果
        
        for idx,translate_word in enumerate(results):
            query_word = word_list[idx]
            mapping[query_word] = translate_word
    return mapping

def load(source_language,target_language,root=DATASET_ROOT):
    from_path = os.path.join(root,"dict_{}_to_{}.pkl".format(source_language,target_language))
    DICT = pd.read_pickle(from_path)
    return DICT


    
if __name__ == "__main__":

    generate("zh","en")