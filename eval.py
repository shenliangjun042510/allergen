import os
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import classification_report, roc_auc_score, roc_curve, precision_recall_curve, auc
from cross_allergen import check_cross_parallel

def evaluate_directory(data_dir, label, use_epitope=False):
    """
    遍历某个文件夹中的FASTA文件，每个文件含两条序列，输出预测结果。
    """
    y_true, y_score = [], []

    fasta_files = sorted([f for f in os.listdir(data_dir) if f.endswith(('.fa', '.fasta'))], key=lambda x: int(x.split('.')[0]))
    print(f"📁 正在处理 {data_dir}，共 {len(fasta_files)} 个样本")

    for file in fasta_files:
        file_path = os.path.join(data_dir, file)
        seqs = list(SeqIO.parse(file_path, "fasta"))

        if len(seqs) != 2:
            print(f"⚠️ 文件 {file} 中不是两个序列，跳过")
            continue

        # 模拟一个小数据库，仅包含 seq2
        db_records = [seqs[1]]

        results = check_cross_parallel(seqs[0], db_records, win=80, id_thresh=0.35, n_jobs=1, local_enabled=use_epitope)
        if results:
            score = results[0]['score']
        else:
            score = 0.0  # 没有匹配上的，设为最低分

        y_true.append(label)
        y_score.append(score)

    return y_true, y_score


def binarize_scores(scores, threshold):
    return [1 if s >= threshold else 0 for s in scores]


def evaluate_all(pos_dir, neg_dir, use_epitope=False, threshold=0.35):
    y_true_pos, y_score_pos = evaluate_directory(pos_dir, label=1, use_epitope=use_epitope)
    y_true_neg, y_score_neg = evaluate_directory(neg_dir, label=0, use_epitope=use_epitope)

    y_true = y_true_pos + y_true_neg
    y_score = y_score_pos + y_score_neg
    y_pred = binarize_scores(y_score, threshold)

    print("\n📊 分类报告：")
    print(classification_report(y_true, y_pred, digits=3))

    auc_score = roc_auc_score(y_true, y_score)
    print(f"🔥 ROC AUC 分数: {auc_score:.3f}")

    # 可视化
    fpr, tpr, _ = roc_curve(y_true, y_score)
    precision, recall, _ = precision_recall_curve(y_true, y_score)

    plt.figure()
    plt.plot(fpr, tpr, label=f'ROC Curve (AUC = {auc_score:.3f})')
    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(f'ROC Curve ({"Epitope" if use_epitope else "WHO"} Based)')
    plt.legend()
    plt.grid(True)
    plt.savefig(f'roc_{"epitope" if use_epitope else "who"}.png')
    print(f"📈 ROC曲线已保存为 roc_{'epitope' if use_epitope else 'who'}.png")

    plt.figure()
    pr_auc = auc(recall, precision)
    plt.plot(recall, precision, label=f'PR Curve (AUC = {pr_auc:.3f})')
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title(f'PR Curve ({"Epitope" if use_epitope else "WHO"} Based)')
    plt.legend()
    plt.grid(True)
    plt.savefig(f'pr_{"epitope" if use_epitope else "who"}.png')
    print(f"📉 PR曲线已保存为 pr_{'epitope' if use_epitope else 'who'}.png")

    return {
        'roc_auc': auc_score,
        'pr_auc': pr_auc
    }


if __name__ == "__main__":
    print("======= 基于 WHO 的模型评估 =======")
    evaluate_all("dataset/positive", "dataset/negative", use_epitope=False, threshold=0.35)

    print("\n======= 基于 Epitope 的模型评估 =======")
    evaluate_all("dataset/positive", "dataset/negative", use_epitope=True, threshold=0.35)
